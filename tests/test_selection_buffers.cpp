#include "tests.h"
#include "bit_info.h"

#include <array>

// Small inputs keep this suitable for a serial, sanitizer-enabled regression.
template <typename T>
void check_selection_buffers() {
    std::array<T, 32> input;
    for (size_t i = 0; i < input.size(); ++i) {
        input[i] = T(1) + T(i % 8) / T(8) + T(i % 3) / T(256);
    }
    const auto original = input;
    std::array<T, 32> output, expected;
    std::string config = test_data_path("real_info.yaml");
    real_info_params<T> params(config.data(), input.data(), input.size());

    for (bool fixed_point : {false, true}) {
        params.fixed_point = fixed_point;
        if (fixed_point) params.fixed_range = get_fixed_point_init_params(input.data(), input.size());
        // Recorded from the existing selectors before replacing their buffers.
        const int expected_bits = fixed_point
            ? (std::is_same_v<T, float> ? 39 : 10)
            : (std::is_same_v<T, float> ? 15 : 44);
        shave_bits(params, input.data(), input.size(), expected_bits, expected.data());

        for (int previous : {0, 1, expected_bits, expected_bits + 2}) {
            for (int repeat = 0; repeat < 3; ++repeat) {
                const int linear = pick_bits_to_shave_template(params, input.data(), input.size(), 0.99, previous);
                test_assert(linear == expected_bits, __FILE__, __LINE__, "Linear selection changed");
                const int binary = pick_bits_to_shave_binary_search_template(
                    params, input.data(), input.size(), 0.99, previous, config.data(), output.data());
                test_assert(binary == expected_bits, __FILE__, __LINE__, "Binary selection changed");
                test_assert(std::memcmp(output.data(), expected.data(), output.size() * sizeof(T)) == 0,
                            __FILE__, __LINE__, "Binary output changed");
                test_assert(preserved_information_template(params, input.data(), output.data(), input.size()) >= 0.99,
                            __FILE__, __LINE__, "Selected output lost information");
            }
        }
        // Exercise the linear selector's fall-through path as well as early returns.
        test_assert(pick_bits_to_shave_template(params, input.data(), input.size(), 1.0, 0) == 0,
                    __FILE__, __LINE__, "Linear fall-through changed");

        std::array<double, 64> redundancy_output{};
        redundancy(params, input.data(), input.data(), input.size(), redundancy_output.data());
        const int bits = fixed_point ? 63 : sizeof(T) * CHAR_BIT;
        for (int bit = 0; bit < bits; ++bit) {
            test_assert(std::abs(redundancy_output[bit] - 1.0) < 1e-12,
                        __FILE__, __LINE__, "Identical arrays lost redundancy");
        }
    }
    test_assert(input == original, __FILE__, __LINE__, "Selection modified input");

    // Invalid data throws after selector scratch buffers have been allocated.
    input[0] = params.fixed_range.maximum + T(1);
    bool linear_rejected = false;
    bool binary_rejected = false;
    try {
        pick_bits_to_shave_template(params, input.data(), input.size(), 0.99, 0);
    } catch (const fixed_point_error &) {
        linear_rejected = true;
    }
    try {
        pick_bits_to_shave_binary_search_template(
            params, input.data(), input.size(), 0.99, 0, config.data(), output.data());
    } catch (const fixed_point_error &) {
        binary_rejected = true;
    }
    test_assert(linear_rejected && binary_rejected, __FILE__, __LINE__, "Invalid range was not rejected");

    // Exercise metric early returns without depending on input-sized scratch.
    std::array<T, 5> zero_information = {T(1), T(1), T(2), T(2), T(1)};
    params.fixed_point = false;
    test_assert(preserved_information_template(params, zero_information.data(), zero_information.data(),
                                              zero_information.size()) == 0.0,
                __FILE__, __LINE__, "Zero-information return changed");
    std::array<T, 3> constant = {T(2), T(2), T(2)};
    std::array<T, 3> constant_output{};
    test_assert(pick_bits_to_shave_template(params, constant.data(), constant.size(), 0.99, 0) == 0,
                __FILE__, __LINE__, "Constant linear selection changed");
    test_assert(pick_bits_to_shave_binary_search_template(params, constant.data(), constant.size(), 0.99, 0,
                                                         config.data(), constant_output.data()) == 0 &&
                constant_output == constant,
                __FILE__, __LINE__, "Constant binary selection changed");
    test_assert(preserved_information_template(params, constant.data(), constant_output.data(), constant.size()) == 1.0,
                __FILE__, __LINE__, "Constant preservation changed");

    std::array<T, 4> patterns = {T(1), T(1), T(2), T(2)};
    test_assert(bitpattern_entropy(patterns.data(), patterns.size()) == 1.0,
                __FILE__, __LINE__, "Entropy sort scratch changed results");
}

int main() {
    check_selection_buffers<float>();
    check_selection_buffers<double>();
    return 0;
}
