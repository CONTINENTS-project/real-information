"""Small loader checks; no native library or NumPy installation is needed."""

import ctypes
import os
from pathlib import Path
import sys
import types
import unittest
from unittest.mock import Mock, patch


SOURCE = Path(__file__).resolve().parents[1] / "src" / "real_info.py"
WRAPPERS = [(SOURCE.read_text(), "lib/bit_info.so")]
if os.environ.get("REAL_INFO_TEST_WRAPPER"):
    WRAPPERS.append((
        Path(os.environ["REAL_INFO_TEST_WRAPPER"]).read_text(),
        os.environ["REAL_INFO_TEST_LIBRARY_RELATIVE_PATH"],
    ))


def import_wrapper(source, prefix, environment, loader, exists=True):
    module = types.ModuleType("real_info_loader_test")
    module.__file__ = str(prefix / "real_info.py")
    with patch.dict(os.environ, environment, clear=True), \
            patch.dict(sys.modules, {"numpy": types.ModuleType("numpy")}), \
            patch.object(Path, "is_file", return_value=exists), \
            patch.object(ctypes, "CDLL", loader):
        exec(compile(source, module.__file__, "exec"), module.__dict__)
    return module


class TestLibraryDiscovery(unittest.TestCase):
    def test_default_follows_wrapper_when_prefix_moves(self):
        for source, relative_library in WRAPPERS:
            for name in ("real-info-install", "real-info-moved install"):
                prefix = Path("/tmp") / name
                with self.subTest(library=relative_library, prefix=name):
                    loader = Mock()
                    module = import_wrapper(source, prefix, {}, loader)
                    expected = (prefix / relative_library).resolve()
                    self.assertEqual(module._lib_path, expected)
                    loader.assert_called_once_with(str(expected))

    def test_explicit_override_wins(self):
        override = Path("/tmp/real-info-override/alternate.so")
        for source, _ in WRAPPERS:
            loader = Mock()
            import_wrapper(source, Path("/tmp/prefix"),
                           {"REAL_INFO_LIB_PATH": str(override)}, loader)
            loader.assert_called_once_with(str(override.resolve()))

    def test_relative_override_is_resolved_before_loading(self):
        for source, _ in WRAPPERS:
            loader = Mock()
            import_wrapper(source, Path("/tmp/prefix"),
                           {"REAL_INFO_LIB_PATH": "alternate.so"}, loader)
            loader.assert_called_once_with(str(Path("alternate.so").resolve()))

    def test_missing_override_does_not_fall_back(self):
        for source, _ in WRAPPERS:
            loader = Mock()
            with self.assertRaisesRegex(FileNotFoundError, "missing.so"):
                import_wrapper(source, Path("/tmp/prefix"),
                               {"REAL_INFO_LIB_PATH": "/tmp/missing.so"},
                               loader, exists=False)
            loader.assert_not_called()

    def test_missing_default_reports_path_and_override(self):
        for source, relative_library in WRAPPERS:
            loader = Mock()
            with self.assertRaises(FileNotFoundError) as caught:
                import_wrapper(source, Path("/tmp/prefix"), {}, loader, exists=False)
            self.assertIn(str((Path("/tmp/prefix") / relative_library).resolve()),
                          str(caught.exception))
            self.assertIn("REAL_INFO_LIB_PATH", str(caught.exception))
            loader.assert_not_called()

    def test_loader_failure_preserves_cause_and_reports_path(self):
        for source, relative_library in WRAPPERS:
            error = OSError("missing OpenMP runtime")
            loader = Mock(side_effect=error)
            with self.assertRaises(OSError) as caught:
                import_wrapper(source, Path("/tmp/prefix"), {}, loader)
            self.assertIn(relative_library, str(caught.exception))
            self.assertIn("missing OpenMP runtime", str(caught.exception))
            self.assertIn("REAL_INFO_LIB_PATH", str(caught.exception))
            self.assertIs(caught.exception.__cause__, error)


if __name__ == "__main__":
    unittest.main()
