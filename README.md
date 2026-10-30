C++ implementation of bitinformation.jl from [Klower et al. 2021](https://www.nature.com/articles/s43588-021-00156-2)

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

The unused `fixed_point_scale_bits` setting has been removed and is rejected in
every scope, including when fixed-point mode is disabled. Fixed-point conversion
continues to map the input range to codes from zero through `2^63 - 1024`; its
precision is not configurable through YAML.

## Bit operations

Shaving, setting, half-shaving, rounding, and grooming accept counts from zero
through 23 for native `float`, 52 for native `double`, or 63 in fixed-point mode.
Native counts modify fraction bits; rounding can carry into the exponent.
Invalid counts throw C++ `fixed_point_error` before output is modified, including
for constant inputs. The C bindings retain their existing fatal-error policy.

All zero-bit array operations copy the input bit patterns exactly, including
signed zeros, without a fixed-point encode/decode round trip. Native operations
preserve infinities and NaN payloads. Finite rounding may overflow to infinity.
Half-shaving clears the lowest `n` bits and sets bit `n - 1`; nearest-even
rounding compares the discarded bits with the midpoint and resolves ties using
the least-significant retained bit. The unsigned implementation backs `bits<T>`
and the `*_bits` array functions. Obsolete transformation and bit-counting
implementations have been removed.

Fixed-point operations keep codes within `[0, Q]`, where `Q = 2^63 - 1024`.
Setting, half-shaving, or rounding above that endpoint saturates to `Q` before
decoding. Constant fixed-point inputs remain unchanged for every valid count.

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
