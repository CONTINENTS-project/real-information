#include "tests.h"
#include "bit_info.h"
#include <array>

template <typename T>
void check_caps() {
    std::array<T, 32> input, output, expected;
    for (size_t i = 0; i < input.size(); ++i)
        input[i] = T(1) + T(i % 8) / T(8) + T(i % 3) / T(256);
    const auto original = input;
    std::string config = test_data_path("real_info.yaml");
    real_info_params<T> params(config.data(), input.data(), input.size());
    for (bool fixed : {false, true}) {
        params.fixed_point = fixed;
        if (fixed) params.fixed_range = get_fixed_point_init_params(input.data(), input.size());
        const int maximum = get_max_shave_bits(params);
        for (int cap : {0, 1, 4, maximum}) {
            for (int previous : {0, 1, maximum + 5}) {
                for (double tolerance : {0.0, 0.99, 1.0}) {
                    output.fill(T(-99));
                    const int selected = pick_bits_to_shave_binary_search_capped_template(
                        params, input.data(), input.size(), tolerance, previous, cap, output.data());
                    test_assert(selected >= 0 && selected <= cap, __FILE__, __LINE__, "Count exceeds cap");
                    shave_bits(params, input.data(), input.size(), selected, expected.data());
                    test_assert(std::memcmp(output.data(), expected.data(), sizeof(output)) == 0,
                                __FILE__, __LINE__, "Selected count and output disagree");
                    if (tolerance == 0.0)
                        test_assert(selected == cap, __FILE__, __LINE__, "Inclusive cap was not tested");
                    if (selected > 0)
                        test_assert(preserved_information_template(params, input.data(), output.data(), input.size()) >= tolerance,
                                    __FILE__, __LINE__, "Selected candidate fails tolerance");
                }
            }
        }
        bool rejected = false;
        try { get_max_shave_bits(params, maximum + 1); }
        catch (const fixed_point_error &) { rejected = true; }
        test_assert(rejected, __FILE__, __LINE__, "Invalid cap accepted");
    }
    test_assert(input == original, __FILE__, __LINE__, "Input was modified");
    input.fill(T(2));
    params.fixed_point = false;
    test_assert(pick_bits_to_shave_binary_search_capped_template(params, input.data(), input.size(),
                0.99, 5, 1, output.data()) == 0 && output == input,
                __FILE__, __LINE__, "Constant output uninitialised");
    for (size_t i = 0; i < input.size(); ++i) input[i] = i % 2 ? T(1e-30) : T(-1e-30);
    params.mask_zero = true;
    test_assert(pick_bits_to_shave_binary_search_capped_template(params, input.data(), input.size(),
                0.99, 5, 1, output.data()) == 0 && output == input,
                __FILE__, __LINE__, "Masked output uninitialised");
}

int main() {
    check_caps<float>();
    check_caps<double>();
    return 0;
}
