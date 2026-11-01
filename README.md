C++ implementation of bitinformation.jl from [Klower et al. 2021](https://www.nature.com/articles/s43588-021-00156-2)

## Build and install

Requirements: CMake 3.16+, a C++20 compiler, and OpenMP 4.5+ (the header uses
array-section reductions). Python usage additionally needs
NumPy. A Fortran compiler is optional; applications compile the installed Fortran
interface source with their own compiler.

Configure from the repository directory with a local installation prefix:

```bash
cmake -S . -B build \
  -DBUILD_TESTING=OFF \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PWD/install"
```

Build and install on a workstation or allocated compute node:

```bash
cmake --build build --parallel 2
cmake --install build
```

The installation contains `include/bit_info.h`, `include/bit_info_mod.f90`,
`real_info.py` at the prefix root, and `lib/bit_info.so`. Standard CMake
install-directory settings apply; for example, `-DCMAKE_INSTALL_LIBDIR=lib64`
puts the library and CMake package under `lib64/`. Keep installation directories
relative to the prefix to allow moving the whole installation.

An external C++ application's `CMakeLists.txt` can use:

```cmake
find_package(real_information CONFIG REQUIRED)
add_executable(my_application main.cpp)
target_link_libraries(my_application PRIVATE real_information::bit_info)
```

Configure that application with `-DCMAKE_PREFIX_PATH=/path/to/install`. The target
supplies the include directory, C++20 requirement, and OpenMP compiler/linker
requirements. CMake package metadata currently uses version `0.1.0`, with version
compatibility within the same minor release. Debugging and optimization options
follow the selected CMake build configuration.

For Python, put the installation prefix on `PYTHONPATH`:

```bash
export PYTHONPATH="$PWD/install${PYTHONPATH:+:$PYTHONPATH}"
python -c 'import real_info'
```

The installed wrapper locates the library relative to itself, including a
configured `lib64/` layout. Moving the prefix preserves that lookup. To use a
source-tree wrapper with a build-tree library, or select another library, set
`REAL_INFO_LIB_PATH` to its full path. 

Fortran applications can compile the installed `include/bit_info_mod.f90` and
link the shared library by its full filename. It has no `lib` prefix, so a plain
`-lbit_info` lookup is insufficient. Arrange the executable's runtime library
search path with the application's CMake target/RPATH settings or, on Linux,
`LD_LIBRARY_PATH` containing the installation's library directory. The installed
CMake target is also available to mixed C++/Fortran CMake applications.


## Configuration

`real_info_params` requires the top-level boolean settings `mantissa_only`,
`mask_zero`, and `fixed_point`. Missing files, required keys, invalid conversions,
and unsupported configuration structures throw C++ `configuration_error`.
Recoverable reporting through the C/Python/Fortran bindings is still to-do after
verification testing; callers must not assume those bindings can recover yet.

The dependency-free reader supports top-level scalar mappings and an optional
`variables` sequence of scalar mappings:

```yaml
mantissa_only: false
mask_zero: true
fixed_point: false

variables:
  - var_name: PRECC
    mantissa_only: true
    mask_zero: false
```

Each entry needs a unique, nonempty `var_name`. Direct C++ `read_yaml_scalar`
calls can select that name; absent overrides and unlisted names inherit only
top-level settings. A null/empty name selects globals. Public bindings do not
accept a variable name yet. Shared GenTS scalar settings such as `real_info`
and `default_real_info` are accepted.

Use spaces for consistent indentation; both indented and indentless sequences
are supported. Use `variables: []` for an empty sequence. Duplicate keys/names,
arbitrary nested mappings/sequences, flow collections, null values, anchors,
tags, and multiline scalars are rejected. This is a defined YAML subset, not a
general-purpose YAML parser.

Strings can be plain, single-quoted (with `''` for an apostrophe), or double-quoted
(with `\"`, `\\`, `\n`, `\r`, and `\t` escapes). A `#` starts a comment outside
quotes when preceded by whitespace or at the start of a scalar. Numeric and
boolean conversions require unquoted values. Booleans accept `true`/`false` and
their title/upper-case spellings. Numbers use decimal notation, optionally
scientific notation for floating point, and must be finite and fit their target
type. Malformed values are rejected rather than partially converted.

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
