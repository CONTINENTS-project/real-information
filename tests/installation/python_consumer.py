"""Exercise a real installed wrapper/library pair from an unrelated directory."""

import os
from pathlib import Path
import sys

import numpy as np
import real_info


prefix, relative_library, config = sys.argv[1:]
expected_library = Path(os.environ.get("REAL_INFO_LIB_PATH",
                                     Path(prefix) / relative_library)).resolve()
assert Path(real_info.__file__).resolve() == Path(prefix) / "real_info.py"
assert real_info._lib_path == expected_library
for dtype in (np.float32, np.float64):
    data = np.full(4, 2.5, dtype=dtype)
    assert real_info.pick_bits_to_shave(data, data.size, 0.01, 0, config) == 0
    count, result = real_info.pick_bits_to_shave_binary_search(
        data, data.size, 0.01, 0, config)
    assert count == 0
    np.testing.assert_array_equal(result, data)
print(f"Installed Python selectors passed: {real_info._lib_path}")
