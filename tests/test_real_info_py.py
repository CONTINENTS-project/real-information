"""
Tests for the input validation of the real_info Python wrapper.

The production entry points ``pick_bits_to_shave`` and
``pick_bits_to_shave_binary_search`` must reject invalid Python-level inputs
(empty arrays, non-float dtypes, non-native byte order, mismatched n_elem)
with normal Python exceptions *before* calling into the C layer, because the
fixed-point C code intentionally terminates the process on invalid input.

This file is a plain unittest module so it can be run either via ctest
(``${Python3_EXECUTABLE} test_real_info_py.py`` as set up in CMakeLists.txt)
or directly with any python that has numpy installed.

Each test is run against every available (real_info.py, bit_info.so)
combination:
  * real-information/src/real_info.py        + build/bit_info.so
  * install/real-information/real_info.py    + install/real-information/lib/bit_info.so
and any cross combinations found.
"""
import importlib.util
import os
import sys
import unittest
from pathlib import Path

try:
    import numpy as np
except ImportError:
    sys.exit("numpy is required to run test_real_info_py.py")

TESTS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = TESTS_DIR.parent
WORKSPACE_ROOT = PROJECT_ROOT.parent

SRC_MODULE = PROJECT_ROOT / "src" / "real_info.py"
INSTALL_MODULE = WORKSPACE_ROOT / "install" / "real-information" / "real_info.py"
CONFIG = TESTS_DIR / "real_info.yaml"  # fixed_point: false


def _lib_candidates(module_path):
    """Ordered, de-duplicated list of usable bit_info.so paths for a module."""
    cands = []
    env = os.environ.get("REAL_INFO_LIB_PATH")
    if env:
        cands.append(Path(env))
    if module_path == SRC_MODULE:
        cands += [PROJECT_ROOT / "build" / "bit_info.so",
                  WORKSPACE_ROOT / "install" / "real-information" / "lib" / "bit_info.so"]
    else:
        cands += [module_path.parent / "lib" / "bit_info.so",
                  PROJECT_ROOT / "build" / "bit_info.so"]
    seen, out = set(), []
    for p in cands:
        p = p.resolve() if p.is_file() else None
        if p is not None and p not in seen:
            seen.add(p)
            out.append(p)
    return out


def _load_module(unique_name, module_path, lib):
    old = os.environ.get("REAL_INFO_LIB_PATH")
    os.environ["REAL_INFO_LIB_PATH"] = str(lib)
    try:
        spec = importlib.util.spec_from_file_location(unique_name, module_path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
    finally:
        if old is None:
            os.environ.pop("REAL_INFO_LIB_PATH", None)
        else:
            os.environ["REAL_INFO_LIB_PATH"] = old
    return mod


def _load_all():
    """Return a list of (label, module) tuples for every usable combination."""
    cases = []
    for i, module_path in enumerate([SRC_MODULE, INSTALL_MODULE]):
        if not module_path.is_file():
            continue
        for j, lib in enumerate(_lib_candidates(module_path)):
            label = f"{module_path.parent.name}/{module_path.name} <- {lib}"
            spec_name = f"real_info_test_target_{i}_{j}"
            cases.append((label, _load_module(spec_name, module_path, lib)))
    return cases


CASES = _load_all()

A32 = (1.0 + 0.12345 * np.arange(16, dtype=np.float32))
A64 = (2.5 + 0.314159 * np.arange(32, dtype=np.float64))
TOL = 0.01


class TestPickSelectorsValidation(unittest.TestCase):
    """Validation behavior of pick_bits_to_shave* across all module builds."""

    def _assert_raises(self, exc_type, label, mod, *args):
        with self.assertRaises(exc_type) as ctx:
            mod.pick_bits_to_shave(*args)
        return ctx

    def setUp(self):
        if not CASES:
            self.skipTest("no (real_info.py, bit_info.so) combination available")

    def test_infer_dtype_raises(self):
        for label, mod in CASES:
            with self.subTest(label):
                self.assertEqual(mod._infer_dtype(A32), "float")
                self.assertEqual(mod._infer_dtype(A64), "double")
                self.assertRaises(ValueError, mod._infer_dtype,
                                  np.array([1, 2, 3]))       # int dtype
                self.assertRaises(ValueError, mod._infer_dtype,
                                  [1.0, 2.0])                # no dtype attr

    def test_valid_float32(self):
        for label, mod in CASES:
            with self.subTest(label):
                bits = mod.pick_bits_to_shave(A32, A32.size, TOL, 0, CONFIG)
                self.assertIsInstance(bits, int)
                self.assertGreaterEqual(bits, 0)

    def test_valid_float64(self):
        for label, mod in CASES:
            with self.subTest(label):
                bits = mod.pick_bits_to_shave(A64, A64.size, TOL, 0, CONFIG)
                self.assertIsInstance(bits, int)
                self.assertGreaterEqual(bits, 0)

    def test_binary_search_returns_bits_and_shaved_array(self):
        for label, mod in CASES:
            with self.subTest(label):
                bits, A_shaved = mod.pick_bits_to_shave_binary_search(
                    A32, A32.size, TOL, 0, CONFIG)
                self.assertIsInstance(bits, int)
                self.assertGreaterEqual(bits, 0)
                self.assertEqual(A_shaved.dtype, np.float32)
                self.assertEqual(A_shaved.shape, A32.shape)

    def test_strided_and_2d_input_accepted(self):
        for label, mod in CASES:
            with self.subTest(label):
                A_strided = A64[::2]  # non-contiguous 1-D view
                self.assertFalse(A_strided.flags["C_CONTIGUOUS"])
                bits = mod.pick_bits_to_shave(A_strided, A_strided.size,
                                              TOL, 0, CONFIG)
                self.assertGreaterEqual(bits, 0)

                A_2d = A32.reshape(2, 8)
                bits = mod.pick_bits_to_shave(A_2d, A_2d.size, TOL, 0, CONFIG)
                self.assertGreaterEqual(bits, 0)

    def test_mismatched_n_elem_raises(self):
        for label, mod in CASES:
            with self.subTest(label):
                self._assert_raises(ValueError, label, mod,
                                    A32, A32.size - 1, TOL, 0, CONFIG)
                self._assert_raises(ValueError, label, mod,
                                    A32, A32.size + 1, TOL, 0, CONFIG)

    def test_empty_array_raises(self):
        for label, mod in CASES:
            with self.subTest(label):
                self._assert_raises(ValueError, label, mod,
                                    np.empty(0, dtype=np.float32), 0,
                                    TOL, 0, CONFIG)

    def test_non_float_dtype_raises(self):
        for label, mod in CASES:
            with self.subTest(label):
                self._assert_raises(ValueError, label, mod,
                                    np.arange(16, dtype=np.int32), 16,
                                    TOL, 0, CONFIG)

    def test_non_native_byte_order_raises(self):
        for label, mod in CASES:
            with self.subTest(label):
                A_be = A32.astype(A32.dtype.newbyteorder(">"))
                self.assertFalse(A_be.dtype.isnative)
                self._assert_raises(TypeError, label, mod,
                                    A_be, A_be.size, TOL, 0, CONFIG)

    def test_misspelled_presereved_information_removed(self):
        for label, mod in CASES:
            with self.subTest(label):
                self.assertFalse(hasattr(mod, "presereved_information"),
                                 "presereved_information should be deleted")


if __name__ == "__main__":
    if not CASES:
        sys.exit("ERROR: no (real_info.py, bit_info.so) combination found; "
                 "build bit_info first (cmake build).")
    print(f"Testing {len(CASES)} real_info module/library combination(s):")
    for label, _ in CASES:
        print(f"  - {label}")
    unittest.main(verbosity=2)
