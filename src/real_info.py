import ctypes
import os
from pathlib import Path
import numpy as np

# Load the compiled C library
_DEFAULT_LIB_PATH = (
    Path(__file__).resolve().parents[2]
    / "install"
    / "real-information"
    / "lib"
    / "bit_info.so"
)

_lib_path = Path(os.environ.get("REAL_INFO_LIB_PATH", _DEFAULT_LIB_PATH))

if not _lib_path.is_file():
    raise FileNotFoundError(
        f"Cannot find bit_info shared library at {_lib_path}. "
        "Set REAL_INFO_LIB_PATH to its full path."
    )

real_info = ctypes.CDLL(str(_lib_path))

def _infer_dtype(A):
    if hasattr(A, "dtype"):
        import numpy as np
        if A.dtype == np.float32:
            return "float"
        elif A.dtype == np.float64:
            return "double"
    return ValueError(f"Unable to infer dtype")

def pick_bits_to_shave(A, n_elem, tol, current_bits, config_path):
    dtype_str = _infer_dtype(A)
    if dtype_str == "float":
        c_func = real_info.pick_bits_to_shave_float
        c_func.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.c_size_t, ctypes.c_double, ctypes.c_int, ctypes.c_char_p]
        c_func.restype = ctypes.c_int
        return c_func(A.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), n_elem, ctypes.c_double(tol), current_bits, config_path.encode('utf-8'))
    elif dtype_str == "double":
        c_func = real_info.pick_bits_to_shave_double
        c_func.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_double, ctypes.c_int, ctypes.c_char_p]
        c_func.restype = ctypes.c_int
        return c_func(A.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n_elem, ctypes.c_double(tol), current_bits, config_path.encode('utf-8'))
    else:
        raise ValueError(f"Unsupported dtype for pick_bits_to_shave: {dtype_str}")

def pick_bits_to_shave_binary_search(A, n_elem, tol, current_bits, config_path):
    dtype_str = _infer_dtype(A)
    if dtype_str == "float":
        c_func = real_info.pick_bits_to_shave_binary_search_float
        c_func.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.c_size_t, ctypes.c_double, ctypes.c_int, ctypes.c_char_p]
        c_func.restype = ctypes.c_int
        return c_func(A.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), n_elem, ctypes.c_double(tol), current_bits, config_path.encode('utf-8'))
    elif dtype_str == "double":
        c_func = real_info.pick_bits_to_shave_binary_search_double
        c_func.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_double, ctypes.c_int, ctypes.c_char_p]
        c_func.restype = ctypes.c_int
        return c_func(A.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n_elem, ctypes.c_double(tol), current_bits, config_path.encode('utf-8'))
    else:
        raise ValueError(f"Unsupported dtype for pick_bits_to_shave_binary_search: {dtype_str}")

def presereved_information(A, B, n_elem, config_path):
    dtype_str = _infer_dtype(A)
    if dtype_str == "float":
        c_func = real_info.preserved_information_float
        c_func.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.c_size_t, ctypes.c_char_p]
        c_func.restype = ctypes.c_double
        return c_func(A.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), B.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), n_elem, config_path.encode('utf-8'))
    elif dtype_str == "double":
        c_func = real_info.preserved_information_double
        c_func.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_char_p]
        c_func.restype = ctypes.c_double
        return c_func(A.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), B.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n_elem, config_path.encode('utf-8'))
    else:
        raise ValueError(f"Unsupported dtype for preserved_information: {dtype_str}")

def shave(A, n_elem, n):
    dtype_str = _infer_dtype(A)
    if dtype_str == "float":
        shaved = np.zeros(n_elem, dtype=np.float32)
        c_func = real_info.shave_float
        c_func.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.c_size_t, ctypes.c_int, ctypes.POINTER(ctypes.c_float)]
        c_func.restype = None
        c_func(A.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), n_elem, n, shaved.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))
        return shaved
    elif dtype_str == "double":
        shaved = np.zeros(n_elem, dtype=np.float64)
        c_func = real_info.shave_double
        c_func.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_int, ctypes.POINTER(ctypes.c_double)]
        c_func.restype = None
        c_func(A.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n_elem, n, shaved.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        return shaved
    else:
        raise ValueError(f"Unsupported dtype for shave: {dtype_str}")

def bit_rounding(A, n_elem, n):
    dtype_str = _infer_dtype(A)
    if dtype_str == "float":
        bit_rounding_result = np.zeros(n_elem, dtype=np.float32)
        c_func = real_info.bit_rounding_float
        c_func.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.c_size_t, ctypes.c_int, ctypes.POINTER(ctypes.c_float)]
        c_func.restype = None
        c_func(A.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), n_elem, n, bit_rounding_result.ctypes.data_as(ctypes.POINTER(ctypes.c_float)))
        return bit_rounding_result
    elif dtype_str == "double":
        bit_rounding_result = np.zeros(n_elem, dtype=np.float64)
        c_func = real_info.bit_rounding_double
        c_func.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_int, ctypes.POINTER(ctypes.c_double)]
        c_func.restype = None
        c_func(A.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n_elem, n, bit_rounding_result.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        return bit_rounding_result
    else:
        raise ValueError(f"Unsupported dtype for bit_rounding: {dtype_str}")
