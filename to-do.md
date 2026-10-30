# Real Information: fixes and cleanup plan

This document preserves the original static code review and the subsequent
discussion of fixed-point conversion. The review inspected only this repository;
it did not build, import, or execute the code. Remaining findings below have not
been revalidated by running tests.

**Verification priority: repair the existing tests (step 12) and establish a
usable validation baseline.** Selected implementation work in steps 2–4 has been
authorized before that baseline; this does not establish that verification is
complete. The numbering below retains the original review order.

- [x] Step 1: fix unsafe fixed-point conversion. Completed in commit `a5f25b8`
  (`Fix unsafe fixed-point conversion and degenerate ranges`).
- [ ] Next: repair the tests described in step 12, including fixture paths,
  representation modes, output buffer sizes, and input-read checks.
- [ ] Continue the remaining correctness and functionality fixes.
- [ ] Follow those repairs with a general tidy and README update.

## 1. Fix unsafe fixed-point conversion — completed

Relevant source: [src/bit_info.h](src/bit_info.h), with error propagation in
[src/bit_info.cpp](src/bit_info.cpp) and [src/real_info.py](src/real_info.py).

The original implementation could round `INT64_MAX` to `2^63` before converting
back to an integer, making the upper endpoint unsafe. Initialization also
dereferenced empty input and divided by zero for constant arrays. Conversion back
could overflow the product of the reciprocal scale and the integer endpoint.

Implemented in `a5f25b8`:

- Map the lower endpoint to zero and the upper endpoint to
  `Q = 2^63 - 1024`, which is exactly representable as both `double` and `int64_t`.
- Store explicit range metadata and use double intermediates for normalization,
  decoding, and differences until the final output conversion.
- Validate finite input values and a finite positive width for nonconstant ranges;
  reject unsupported ranges explicitly.
- Reject empty initialization and analysis inputs. Empty elementwise shaving is a
  no-op.
- Handle constants without scaling: code zero decodes to the stored constant,
  fixed-point array operations preserve the input, and automatic selection returns
  zero bits shaved. Copying the input preserves signed-zero patterns.
- Use an explicit constant-reference preservation convention: `1.0` for identical
  input/output bit patterns and `0.0` otherwise.
- Construct encoded zero directly, avoiding an attempt to encode floating-point
  zero when zero lies outside the initialized range.
- Keep fixed-point validation exceptions inside the C++ boundary. C callers can
  inspect `real_info_last_fixed_point_error()`; Python translates these errors to
  `ValueError`. Existing C/Python call signatures remain unchanged.

C++ fixed-point initialization and constructors now take explicit range metadata.
The fixed-point arithmetic operators were adapted to that construction, addressing
the original review's nonexistent-default-constructor finding.

Regression source was added in
[tests/test_fixed_point_conversion.cpp](tests/test_fixed_point_conversion.cpp)
for endpoint mapping, small and wide ranges, constants, signed zeros, empty and
nonfinite input, and C error reporting. It was not built or executed by the
assistant. Step 1 is complete as an implementation change; verification remains
part of repairing and establishing the test baseline.

## 2. Repair the Fortran ABI — implemented

Relevant source: [src/bit_info_mod.f90](src/bit_info_mod.f90) and
[src/bit_info.cpp](src/bit_info.cpp).

Fortran declarations omit configuration arguments required by C++, and
binary-search declarations also omit the output array. Lengths use `c_int` where
C++ expects `size_t`.

- [x] Match every declaration to its implementation.
- [x] Use `c_size_t` for array lengths.
- [x] Provide interoperable configuration strings and output buffers.
- [x] Resolve the unsupplied `kinds`/`confirm` dependencies.
- [ ] Run `test_fortran_abi` against a newly built library in an approved
  compute allocation. The Fortran sources compile, but the available existing
  shared library predates the current source.

## 3. Validate Python arrays before passing pointers — implemented

Relevant source: [src/real_info.py](src/real_info.py) (mirrored in
`install/real-information/real_info.py`).

Scope refinement: validation is only needed for the two production selectors,
`pick_bits_to_shave` and `pick_bits_to_shave_binary_search` (the latter is also
used by GenTS). All other wrappers exist for testing purposes and are not part
of the production data path, so they keep their existing behavior.

- [x] Add a shared `_validated_input(A, n_elem)` helper used by both
  `pick_bits_to_shave*` wrappers. It converts array-likes to numpy, rejects
  non-float dtypes (`TypeError`/`ValueError` from `_infer_dtype`), rejects
  empty arrays and non-native byte order, converts strided (non-contiguous)
  arrays to C-contiguous copies, and rejects mismatched `n_elem`. This also
  protects against the C layer terminating the process on invalid input.
- [x] Route both `pick_bits_to_shave` and `pick_bits_to_shave_binary_search`
  through the helper (both precisions, both the src and the installed copy).
- [x] Make `_infer_dtype` raise `ValueError` instead of returning one.
- [x] Remove the misspelled `presereved_information` wrapper from both copies;
  no call sites existed, so no other code needed corrections.
- [x] Add [tests/test_real_info_py.py](tests/test_real_info_py.py), which is
  auto-registered by the CMake test glob and exercises valid inputs,
  strided/2-D inputs, and each rejection path against every available
  (real_info.py, bit_info.so) combination.

Note: CMake runs Python tests with its discovered `Python3_EXECUTABLE`; if that
interpreter lacks numpy the test is registered but fails, so configure with
`-DPython3_EXECUTABLE=<python-with-numpy>` (e.g.
`/glade/work/aborissov/env/bin/python`).

## 4. Make configuration failures explicit

Relevant source: `read_yaml_scalar` and `real_info_params` in
[src/bit_info.h](src/bit_info.h).

Parameter constructors now initialize their fields and load the configuration
source once per construction. Missing files/required keys and invalid scalar
conversions throw `configuration_error` with source context. Numeric conversions
consume the entire value, enforce the destination range, and reject nonfinite
floating-point values. These are implementation changes; full verification is
still pending.

- [x] Initialize fields and validate conversions; add focused regression source
  in [tests/test_config_validation.cpp](tests/test_config_validation.cpp).
  This test and `test_read_yaml_scalar` were compiled individually and passed;
  the complete build/CTest validation baseline remains pending.
- [ ] TODO after all verification testing is complete: report configuration
  errors through the C, Python, and Fortran bindings. Until then, do not rely on
  recoverable configuration errors through those bindings.
- [ ] Enforce global and variable configuration scopes.
- [ ] TODO after all verification testing is complete: connect variable overrides
  to the public API. Existing exported entry points still use global settings.
- [ ] Implement or remove the unused `fixed_point_scale_bits` setting.

## 5. Remove leaks and input-sized stack allocations

Relevant source: `redundancy`, `preserved_information_template`, and the selection
functions in [src/bit_info.h](src/bit_info.h).

Redundancy and preservation calculations leak their allocated arrays. Linear
selection also leaks its scratch array on successful early returns. Recursive
binary search allocates two input-sized variable-length arrays per call; these
arrays also conflict with the project's standard C++20 configuration.

- [ ] Replace manual allocations with owned, reusable buffers.
- [ ] Release resources on every return path.
- [ ] Replace recursive binary search with an iterative implementation.

## 6. Correct rounding and half-shaving

Relevant source: rounding and half-shaving helpers in
[src/bit_info.h](src/bit_info.h).

Bit rounding uses a discarded bit for tie handling rather than the last retained
bit, and rounding away zero bits performs a negative shift. Half-shaving uses
incorrect bit positions and a 32-bit starting mask for doubles.

- [ ] Consolidate the operations around unsigned integer representations.
- [ ] Make zero-bit operations identities.
- [ ] Validate bit counts.
- [ ] Implement consistent nearest-even rounding.

## 7. Use consistent samples for information metrics

Relevant source: pair counting, marginal entropy, and redundancy in
[src/bit_info.h](src/bit_info.h).

Pair counting excludes a pair only when both values are masked, whereas marginal
entropy excludes masked values independently. Normalized redundancy consequently
combines different populations.

- [ ] Use one eligibility mask for each paired calculation.
- [ ] Replace `M[i] == HA[i] == HB[i] == 0` with explicit conditions.
- [ ] Initialize every output.
- [ ] Define remaining no-sample and zero-information cases consistently with the
  constant-array convention introduced in step 1.

## 8. Repair selection boundaries and output guarantees

Relevant source: the linear and binary selectors in
[src/bit_info.h](src/bit_info.h).

The selectors disagree on strict versus inclusive tolerance. Fresh binary searches
never evaluate the maximum permitted count; linear searches resumed from a
previous count can enter exponent bits. The general all-masked binary-search
return leaves the output buffer untouched. Step 1 addressed the constant-array
case, including filling its output, but the remaining cases still need repair.

- [ ] Use common inclusive bounds.
- [ ] Validate tolerance and previous bit counts.
- [ ] Ensure the returned count always describes the returned array.
- [ ] Establish monotonicity before relying on binary search for fixed-point and
  masking modes.

## 9. Handle entropy edge cases correctly

Relevant source: `bitpattern_entropy` and `bitwise_real_information` in
[src/bit_info.h](src/bit_info.h).

Bit-pattern entropy sorts floating-point values rather than their bit patterns,
merging signed zeros and mishandling NaNs. Empty input also underflows
`n_elem - 1`; adjacent information has the same problem. Step 1's initialization
and selector guards do not replace checks in these individual routines.

- [ ] Count unsigned bit representations.
- [ ] Define behavior for arrays with fewer than two elements.

## 10. Repair build and installation support

Relevant source: [CMakeLists.txt](CMakeLists.txt),
[src/bit_info.h](src/bit_info.h), and [src/real_info.py](src/real_info.py).

CMake hardcodes OpenMP flags without propagating them to consumers of
header-defined algorithms. The core header requires an apparently unused `mpi.h`,
but CMake declares no MPI dependency. Installation omits public headers, and Python
library discovery assumes a surrounding checkout layout.

- [ ] Remove the unused MPI dependency.
- [ ] Use CMake's OpenMP target.
- [ ] Install the public interfaces.
- [ ] Locate the shared library relative to its installation.

## 11. Bring examples into agreement with the API

Relevant source: [examples/scaling.py](examples/scaling.py),
[examples/fixed_point.py](examples/fixed_point.py), and the C++ examples.

Python scaling omits required configuration arguments, treats a tuple result as an
integer, and calls `preserved_information`, while the wrapper exposes the
misspelled `presereved_information`. Fixed-point analysis additionally calls
missing information/frequency bindings, swaps result fields when writing rows,
and assumes matching floating/fixed bit widths. C++ examples use outdated
constructor/helper signatures and must also account for step 1's range metadata.

- [ ] Repair example calls and their return-value handling.
- [ ] Repair TSV schemas, result-field ordering, and plot coordinates.
- [ ] Update descriptions to match the actual calculations.

## 12. Repair test sources and plan missing coverage — do this next

Relevant source: [tests/CMakeLists.txt](tests/CMakeLists.txt),
[tests/test_yaml_override.yaml](tests/test_yaml_override.yaml),
[tests/tests.h](tests/tests.h), and individual test sources.

The shared configuration enables 63-bit fixed-point output, while several tests
allocate only 32 slots and compare floating-point expectations. Tests also use
machine-specific absolute paths, and input loading does not check read failures.

- [ ] Separate floating-point and fixed-point fixtures and expected results.
- [ ] Size output buffers for the selected representation and bit range.
- [ ] Replace machine-specific paths with portable fixture paths supplied through
  CMake or derived from the test source location.
- [ ] Check fixture files open successfully and contain the expected input data.
- [ ] Include the new step 1 regression source in the validation baseline.
- [ ] Plan coverage for boundaries, rounding, configuration, ABI, and array layout.
- [ ] Once execution is authorized, establish a usable baseline before further
  implementation work. Heavy validation on NCAR systems requires an approved
  compute allocation.

## General tidy after correctness and functionality repairs

- [ ] Consolidate overlapping rounding implementations.
- [ ] Separate configuration, metrics, selection, and bindings.
- [ ] Introduce consistent namespaces and immutable constants.
- [ ] Remove debug output and obsolete comments.
- [ ] Address duplicate definitions in `src/my_rounding.h`.
- [ ] Address the non-inline header definition of `entropy`.
- [ ] Replace unsafe `std::bitset` reinterpretation.
- [ ] Keep generated binaries and caches in build/output directories.

The original finding about fixed-point operators attempting to use a nonexistent
default constructor was addressed while adapting conversion in step 1.

Before finalizing documentation, decide whether preservation should continue to
force floating-point analysis when fixed-point mode is selected. Document that
choice or change it, and make example threshold plots use the same definition.

## README update after the build repairs

- [ ] Replace the current brief [README.md](README.md) with an informative
  repository description, installation instructions, and usage documentation.

The following is proposed README material. The commands are for future
documentation; they have not been executed as part of this plan.

### Real Information

A C++20 library for analyzing information in floating-point data and selecting
precision reduction through bit shaving and rounding. It computes per-bit entropy,
mutual information between neighboring samples, redundancy, and estimated
information preservation. It supports floating-point analysis and an experimental
fixed-point representation.

Repository contents:

- `src/`: C++ algorithms, C entry points, Python bindings, Fortran interface
  declarations, rounding utilities, and CLI source.
- `examples/`: synthetic-data scaling and fixed-point comparisons, YAML
  configurations, and example results.
- `tests/`: C++ test sources and reference fixtures registered through CTest.

Requirements: CMake 3.16+, a C++20 compiler with OpenMP support, and Python with
NumPy for the Python interface. Plotting examples additionally require Matplotlib
and Python 3.9+.

From the repository directory, install into a local prefix. Run the build on a
workstation or allocated compute node, as appropriate:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install numpy

cmake -S . -B build \
  -DBUILD_TESTING=OFF \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PWD/install"

cmake --build build --parallel 2
cmake --install build

export PYTHONPATH="$PWD/install${PYTHONPATH:+:$PYTHONPATH}"
export REAL_INFO_LIB_PATH="$PWD/install/lib/bit_info.so"
```

On NSF NCAR systems, load suitable compiler/CMake modules and perform the build on
an approved compute allocation. Obtain approval before submitting this replacement
for the build command:

```bash
qcmd -A NTDD0005 -- cmake --build build --parallel 2
```

Configuration documentation should explain `mantissa_only`, `mask_zero`, and
`fixed_point`; bit ordering and output lengths; preservation semantics; and
supported input shapes and dtypes. Include a minimal Python example once the
bindings are repaired, plus the status and build instructions for the Fortran
interface and CLI.
