C++ implementation of bitinformation.jl from [Klower et al. 2021](https://www.nature.com/articles/s43588-021-00156-2)

## Fortran interface

Compile `src/bit_info_mod.f90` with the same Fortran compiler used for the
application and link against `bit_info.so`. CMake also installs the module
source in `include/`. Array lengths must have kind `c_size_t`, bit counts must
have kind `c_int`, and configuration paths must be NUL-terminated C character
arrays. For example:

```fortran
program example
	use iso_c_binding, only: c_float, c_int, c_size_t, c_char, c_null_char
	use bit_info, only: shave
	implicit none
	character(len=*), parameter :: path = 'config.yaml'
	character(kind=c_char) :: config(len(path) + 1)
	real(c_float) :: data(3), result(3)
	integer :: index

	config = [(path(index:index), index = 1, len(path)), c_null_char]
	data = [1.0_c_float, 1.5_c_float, 2.0_c_float]
	call shave(data, size(data, kind=c_size_t), 2_c_int, result, config)
end program
```

The binary-search selector additionally takes an output array after the
configuration path. Output arrays must have at least the requested number of
elements. Fixed-point validation failures terminate the calling process.
When testing is enabled, CMake registers `test_fortran_abi` if it finds a
Fortran compiler.
