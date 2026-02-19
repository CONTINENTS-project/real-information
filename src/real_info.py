import ctypes

def _infer_dtype(A):
    if hasattr(A, "dtype"):
        import numpy as np
        if A.dtype == np.float32:
            return "float"
        elif A.dtype == np.float64:
            return "double"
    return ValueError(f"Unable to infer dtype")

def pick_bits_to_shave(A, n_elem, tol, current_bits):
    dtype_str = _infer_dtype(A)
    if dtype_str == "float":
        c_func = real_info.pick_bits_to_shave_float
        c_func.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.c_size_t, ctypes.c_double, ctypes.c_int]
        c_func.restype = ctypes.c_int
        return c_func(A.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), n_elem, ctypes.c_double(tol), current_bits)
    elif dtype_str == "double":
        c_func = real_info.pick_bits_to_shave_double
        c_func.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_double, ctypes.c_int]
        c_func.restype = ctypes.c_int
        return c_func(A.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n_elem, ctypes.c_double(tol), current_bits)
    else:
        raise ValueError(f"Unsupported dtype for pick_bits_to_shave: {dtype_str}")

def pick_bits_to_shave_binary_search(A, n_elem, tol, current_bits):
    dtype_str = _infer_dtype(A)
    if dtype_str == "float":
        c_func = real_info.pick_bits_to_shave_binary_search_float
        c_func.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.c_size_t, ctypes.c_double, ctypes.c_int]
        c_func.restype = ctypes.c_int
        return c_func(A.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), n_elem, ctypes.c_double(tol), current_bits)
    elif dtype_str == "double":
        c_func = real_info.pick_bits_to_shave_binary_search_double
        c_func.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_double, ctypes.c_int]
        c_func.restype = ctypes.c_int
        return c_func(A.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n_elem, ctypes.c_double(tol), current_bits)
    else:
        raise ValueError(f"Unsupported dtype for pick_bits_to_shave_binary_search: {dtype_str}")

def presereved_information(A, B, n_elem):
    dtype_str = _infer_dtype(A)
    if dtype_str == "float":
        c_func = real_info.pick_bits_to_shave_binary_search_float
        c_func.argtypes = [ctypes.POINTER(ctypes.c_float), ctypes.POINTER(ctypes.c_float), ctypes.c_size_t]
        c_func.restype = ctypes.c_double
        return c_func(A.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), B.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), n_elem)
    elif dtype_str == "double":
        c_func = real_info.pick_bits_to_shave_binary_search_double
        c_func.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_size_t]
        c_func.restype = ctypes.c_double
        return c_func(A.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), B.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), n_elem)
    else:
        raise ValueError(f"Unsupported dtype for preserved_information: {dtype_str}")

def shave(A, n_elem, n, shaved):
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
