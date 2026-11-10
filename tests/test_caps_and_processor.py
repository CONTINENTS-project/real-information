"""Small serial cap/timing checks; all data are synthetic and fit in kilobytes."""
import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]


class CapsAndTimingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        library = Path(os.environ.get('REAL_INFO_LIB_PATH', ROOT / 'real-information/build/bit_info.so'))
        if not library.is_file():
            raise unittest.SkipTest('Build bit_info first, or set REAL_INFO_LIB_PATH')
        with patch.dict(os.environ, REAL_INFO_LIB_PATH=str(library.resolve())):
            spec = importlib.util.spec_from_file_location('real_info', ROOT / 'real-information/src/real_info.py')
            cls.native = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(cls.native)
        cls.previous_module = sys.modules.get('real_info')
        sys.modules['real_info'] = cls.native
        sys.path.insert(0, str(ROOT / 'GenTS'))
        from gents.real_info_processor import RealInfoProcessor
        cls.Processor = RealInfoProcessor

    @classmethod
    def tearDownClass(cls):
        if cls.previous_module is not None:
            sys.modules['real_info'] = cls.previous_module
        else:
            sys.modules.pop('real_info', None)

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.config = self.root / 'config.yaml'
        self.options = dict(use_real_info=True, default_real_info=0.0,
            real_info_eval_freq=2, use_natural_ordering=False,
            mantissa_only=False, mask_zero=False, fixed_point=False,
            max_shave_bits=8, variables=[{'var_name': 'U', 'max_shave_bits': 1}])
        values = 1 + np.arange(32, dtype=np.float32) / 127
        self.data = np.tile(values, (3, 2, 1))
        data = self.data

        class Dataset:
            def get_var_dtype(self, variable): return data.dtype
            def get_var_data_shape(self, variable): return list(data.shape)
            def get_var_dimensions(self, variable): return ['time', 'lev', 'ncol']
            def get_var_attrs(self, variable): return {'units': 'K'}
            def get_global_attrs(self): return {'title': 'synthetic cap test'}
            def get_var_vals(self, variable, time_index_start=None, time_index_end=None):
                if variable == 'lat': return np.linspace(-80, 80, 32)
                if variable == 'lon': return np.arange(32) * 360 / 32
                return data if time_index_start is None else data[time_index_start:time_index_end]
        self.dataset = Dataset()

    def processor(self):
        self.config.write_text(yaml.safe_dump(self.options))
        return self.Processor(self.config)

    def test_native_wrapper_bounds_and_validation(self):
        self.processor()
        for dtype, maximum in ((np.float32, 23), (np.float64, 52)):
            data = self.data[0, 0].astype(dtype)
            for cap in (0, 1, 4, maximum):
                bits, output = self.native.pick_bits_to_shave_binary_search(
                    data, data.size, 0.0, 100, self.config, max_shave_bits=cap)
                self.assertEqual(bits, cap)
                expected = self.native.shave(data, data.size, cap, str(self.config))
                np.testing.assert_array_equal(output.view('u' + str(data.dtype.itemsize)),
                                              expected.view('u' + str(data.dtype.itemsize)))
            with self.assertRaises(ValueError):
                self.native.pick_bits_to_shave_binary_search(data, data.size, .99, 0,
                    self.config, max_shave_bits=maximum + 1)
        for invalid in (True, -1, 1.5, 64):
            with self.assertRaises(ValueError):
                self.native.pick_bits_to_shave_binary_search(self.data[0, 0], 32, .99, 0,
                    self.config, max_shave_bits=invalid)

    def test_per_variable_caps_cached_counts_and_whole_call_timing(self):
        processor = self.processor()
        previous = np.full(2, 20, dtype=np.int32)
        with patch.object(self.native, 'pick_bits_to_shave_binary_search', wraps=self.native.pick_bits_to_shave_binary_search) as search, \
                patch.object(self.native, 'shave', wraps=self.native.shave) as shave:
            _, counts = processor.shave_data(self.data[:1], self.dataset, 'U',
                ['time', 'lev', 'ncol'], previous, 0)
            np.testing.assert_array_equal(counts, [1, 1])
            np.testing.assert_array_equal(previous, [20, 20])
            self.assertEqual(search.call_count, 2)
            self.assertTrue(all(call.kwargs['max_shave_bits'] == 1 for call in search.call_args_list))
            shave.assert_not_called()  # Search returns the capped output directly.
            processor.shave_data(self.data[:1], self.dataset, 'U',
                ['time', 'lev', 'ncol'], previous, 1)
            self.assertEqual(shave.call_count, 2)
            self.assertTrue(all(call.args[2] == 1 for call in shave.call_args_list))
        summary = processor.timing_summary('U')
        self.assertEqual(summary['modes']['selection']['calls'], 1)
        self.assertEqual(summary['modes']['cached']['calls'], 1)
        self.assertGreaterEqual(summary['seconds'], 0)
        _, counts = processor.shave_data(self.data[:1], self.dataset, 'T',
            ['time', 'lev', 'ncol'], np.zeros(2, dtype=np.int32), 0)
        np.testing.assert_array_equal(counts, [8, 8])

    def test_fixed_counts_zero_caps_and_disabled_mode(self):
        self.options.update(shave_bits=63, max_shave_bits=0)
        processor = self.processor()
        for variable, expected in (('U', 1), ('T', 0)):
            output, bits = processor.shave_data(self.data[:1], self.dataset, variable,
                ['time', 'lev', 'ncol'], np.zeros(2, dtype=np.int32), 0)
            np.testing.assert_array_equal(bits, [expected, expected])
            self.assertEqual(processor.timing_summary(variable)['modes']['fixed']['calls'], 1)
            if expected == 0: np.testing.assert_array_equal(output, self.data[:1])
        disabled = self.Processor()
        output, _ = disabled.shave_data(self.data[:1], self.dataset, 'T',
            ['time', 'lev', 'ncol'], [0, 0], 0)
        np.testing.assert_array_equal(output, self.data[:1])
        self.assertEqual(disabled.last_timing['mode'], 'bypass')

    def test_cap_validation_and_failed_timing(self):
        for invalid in (True, -1, 1.5, 64):
            self.options['max_shave_bits'] = invalid
            with self.assertRaises(ValueError): self.processor()
        self.options['max_shave_bits'] = 24
        processor = self.processor()
        with self.assertRaisesRegex(ValueError, 'representation limit'):
            processor.shave_data(self.data[:1], self.dataset, 'T',
                ['time', 'lev', 'ncol'], [0, 0], 0)
        self.assertTrue(processor.last_timing['failed'])
        self.assertEqual(processor.timing_summary('T')['modes']['selection']['failed_calls'], 1)

    def test_output_metadata_counts_and_worker_summary(self):
        from gents.timeseries import write_timeseries_file, generate_time_series
        import netCDF4
        processor = self.processor()
        output = self.root / 'U.nc'
        write_timeseries_file(self.dataset, str(output), 'U', {}, real_info_processor=processor)
        with netCDF4.Dataset(output) as dataset:
            summary = json.loads(dataset['U'].gents_shave_timing)
            self.assertEqual(summary['calls'], 3)
            self.assertEqual(summary['modes']['selection']['calls'], 2)
            self.assertEqual(summary['modes']['cached']['calls'], 1)
            self.assertEqual(dataset['U'].max_shave_bits, 1)
            np.testing.assert_array_equal(dataset['bits_shaved'][:], np.ones((3, 2)))
        # Reusing a file must not falsely attribute old computation to this run.
        write_timeseries_file(self.dataset, str(output), 'U', {}, real_info_processor=processor)
        self.assertEqual(processor.timing_summary('U'), {})
        from contextlib import nullcontext
        with patch('gents.timeseries.MHFDataset', return_value=nullcontext(self.dataset)):
            result = generate_time_series([], str(self.root / 'case'), {},
                {'U': {'ts_string': 'period'}}, processor, return_timings=True)
        self.assertEqual(len(result['paths']), 1)
        self.assertEqual(result['timings'][0]['status'], 'computed')
        self.assertEqual(result['timings'][0]['calls'], 3)


if __name__ == '__main__':
    unittest.main()
