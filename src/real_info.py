import ctypes
import os
from pathlib import Path
import numpy as np

# CMake substitutes the relative path in installed copies. Unconfigured source
# copies use the conventional layout; set REAL_INFO_LIB_PATH for a build-tree .so.
_LIBRARY_RELATIVE_PATH = "@REAL_INFO_LIBRARY_RELATIVE_PATH@"
if _LIBRARY_RELATIVE_PATH.startswith("@"):
    _LIBRARY_RELATIVE_PATH = "lib/bit_info.so"
_DEFAULT_LIB_PATH = Path(__file__).resolve().parent / _LIBRARY_RELATIVE_PATH

_lib_path = Path(os.environ.get("REAL_INFO_LIB_PATH", _DEFAULT_LIB_PATH)).expanduser().resolve()

if not _lib_path.is_file():
    raise FileNotFoundError(
        f"Cannot find bit_info shared library at {_lib_path}. "
        "Set REAL_INFO_LIB_PATH to its full path."
    )

try:
    real_info = ctypes.CDLL(str(_lib_path))
except OSError as error:
    raise OSError(
        f"Cannot load bit_info shared library at {_lib_path}: {error}. "
        "Check its runtime dependencies or set REAL_INFO_LIB_PATH to the full "
        "path of a compatible library."
    ) from error

real_info.real_info_last_fixed_point_error.argtypes = []
real_info.real_info_last_fixed_point_error.restype = ctypes.c_char_p


def _check_fixed_point_result(result, function, arguments):
    error = real_info.real_info_last_fixed_point_error()
    if error is not None:
        raise ValueError(error.decode("utf-8"))
    return result


# Check only the entry points that initialize or validate fixed-point ranges.
# The existing Python return values remain unchanged.
for _function_name in (
    "shave_float", "shave_double",
    "halfshave_float", "halfshave_double",
    "bit_rounding_float", "bit_rounding_double",
    "set_float", "set_double",
    "groom_float", "groom_double",
    "preserved_information_float", "preserved_information_double",
    "pick_bits_to_shave_float", "pick_bits_to_shave_double",
    "pick_bits_to_shave_binary_search_float", "pick_bits_to_shave_binary_search_double",
):
    getattr(real_info, _function_name).errcheck = _check_fixed_point_result
del _function_name

def _infer_dtype(A):
    if hasattr(A, "dtype"):
        if A.dtype == np.float32:
            return "float"
        elif A.dtype == np.float64:
            return "double"
    raise ValueError(f"Unable to infer dtype for {type(A).__name__}")

def _validated_input(A, n_elem):
    """
    Validate input for the production ``pick_bits_to_shave*`` wrappers, so
    that invalid Python-level inputs raise a normal Python exception
    instead of failing inside the C layer (the fixed-point C code
    intentionally terminates the process on invalid input).

    Checks performed:
      * the input is (or can be converted to) a numpy array of float32 or
        float64 values,
      * the array is non-empty,
      * the dtype has native byte order,
      * ``n_elem`` matches the total number of elements.

    A non-contiguous array is accepted and converted to a C-contiguous
    array internally.

    :param A: input array or array-like.
    :param n_elem: number of elements the caller claims the array holds.
    :return: tuple ``(A, dtype_str, n_elem)`` where ``A`` is guaranteed to
             be a C-contiguous float32/float64 numpy array and
             ``dtype_str`` is either ``"float"`` or ``"double"``.
    :raises ValueError: if the dtype is not float32/float64, the array is
            empty, or ``n_elem`` does not equal ``A.size``.
    :raises TypeError: if the array byte order is not native.
    """
    if not isinstance(A, np.ndarray):
        A = np.array(A)
    if not A.dtype.isnative:
        raise TypeError(
            f"Input array has non-native byte order ({A.dtype}); "
            "convert it first, e.g. A = A.astype(A.dtype.newbyteorder('='))"
        )
    dtype_str = _infer_dtype(A)  # raises ValueError for non-float32/64
    if A.size == 0:
        raise ValueError("Input array must be non-empty")
    if not A.flags["C_CONTIGUOUS"]:
        A = np.ascontiguousarray(A)
    if int(n_elem) != A.size:
        raise ValueError(f"n_elem ({int(n_elem)}) does not match A.size ({A.size})")
    return A, dtype_str, int(n_elem)

def pick_bits_to_shave(A, n_elem, tol, current_bits, config_path):
    """
    Determine how many additional bits can be shaved from the mantissa of
    A while keeping the relative information loss within ``tol``.

    :param A: array of float32 or float64 values.
    :param n_elem: number of elements in A (must equal ``A.size``).
    :param tol: tolerable relative information loss.
    :param current_bits: number of bits already shaved from A.
    :param config_path: path to yaml config file.
    :return: number of bits that can still be shaved (int).
    """
    if isinstance(config_path, Path):
        config_path = str(config_path)
    A, dtype_str, n_elem = _validated_input(A, n_elem)
    if dtype_str == "float":
        c_func = real_info.pick_bits_to_shave_float
        c_func.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.c_size_t, ctypes.c_double, ctypes.c_int, ctypes.c_char_p]
        c_func.restype = ctypes.c_int
        return c_func(A.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), n_elem, ctypes.c_double(tol), current_bits, config_path.encode('utf-8'))
    else:
        c_func = real_info.pick_bits_to_shave_double
        c_func.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_double, ctypes.c_int, ctypes.c_char_p]
        c_func.restype = ctypes.c_int
        return c_func(A.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n_elem, ctypes.c_double(tol), current_bits, config_path.encode('utf-8'))

def pick_bits_to_shave_binary_search(A, n_elem, tol, current_bits, config_path):
    """
    Determine how many additional bits can be shaved from the mantissa of
    A using binary search, and return the shaved copy of A.

    :param A: array of float32 or float64 values.
    :param n_elem: number of elements in A (must equal ``A.size``).
    :param tol: tolerable relative information loss.
    :param current_bits: number of bits already shaved from A.
    :param config_path: path to yaml config file.
    :return: tuple ``(n_bits, A_shaved)`` where ``A_shaved`` is a contiguous
             numpy array of the same dtype as the input.
    """
    if isinstance(config_path, Path):
        config_path = str(config_path)
    A, dtype_str, n_elem = _validated_input(A, n_elem)
    if dtype_str == "float":
        c_func = real_info.pick_bits_to_shave_binary_search_float
        c_func.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.c_size_t, ctypes.c_double, ctypes.c_int, ctypes.c_char_p, ctypes.POINTER(ctypes.c_float)]
        c_func.restype = ctypes.c_int
        A_shaved = np.zeros(n_elem, dtype=np.float32)
        return (c_func(A.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), n_elem, ctypes.c_double(tol), current_bits, config_path.encode('utf-8'), A_shaved.ctypes.data_as(ctypes.POINTER(ctypes.c_float))), A_shaved)
    else:
        c_func = real_info.pick_bits_to_shave_binary_search_double
        c_func.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_double, ctypes.c_int, ctypes.c_char_p, ctypes.POINTER(ctypes.c_double)]
        c_func.restype = ctypes.c_int
        A_shaved = np.zeros(n_elem, dtype=np.float64)
        return (c_func(A.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n_elem, ctypes.c_double(tol), current_bits, config_path.encode('utf-8'), A_shaved.ctypes.data_as(ctypes.POINTER(ctypes.c_double))), A_shaved)

def shave(A, n_elem, n, config_path):
    dtype_str = _infer_dtype(A)
    if dtype_str == "float":
        shaved = np.zeros(n_elem, dtype=np.float32)
        c_func = real_info.shave_float
        c_func.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.c_size_t, ctypes.c_int, ctypes.POINTER(ctypes.c_float), ctypes.c_char_p]
        c_func.restype = None
        c_func(A.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), n_elem, n, shaved.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), config_path.encode('utf-8'))
        return shaved
    elif dtype_str == "double":
        shaved = np.zeros(n_elem, dtype=np.float64)
        c_func = real_info.shave_double
        c_func.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.c_char_p]
        c_func.restype = None
        c_func(A.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n_elem, n, shaved.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), config_path.encode('utf-8'))
        return shaved
    else:
        raise ValueError(f"Unsupported dtype for shave: {dtype_str}")

def bit_rounding(A, n_elem, n, config_path):
    dtype_str = _infer_dtype(A)
    if dtype_str == "float":
        bit_rounding_result = np.zeros(n_elem, dtype=np.float32)
        c_func = real_info.bit_rounding_float
        c_func.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.c_size_t, ctypes.c_int, ctypes.POINTER(ctypes.c_float), ctypes.c_char_p]
        c_func.restype = None
        c_func(A.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), n_elem, n, bit_rounding_result.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), config_path.encode('utf-8'))
        return bit_rounding_result
    elif dtype_str == "double":
        bit_rounding_result = np.zeros(n_elem, dtype=np.float64)
        c_func = real_info.bit_rounding_double
        c_func.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.c_char_p]
        c_func.restype = None
        c_func(A.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n_elem, n, bit_rounding_result.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), config_path.encode('utf-8'))
        return bit_rounding_result
    else:
        raise ValueError(f"Unsupported dtype for bit_rounding: {dtype_str}")
