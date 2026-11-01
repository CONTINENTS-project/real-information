#include "tests.h"
#include "bit_info.h"

using MetricOutput = std::array<double, 65>;

MetricOutput unwritten_output() {
    MetricOutput output;
    output.fill(std::numeric_limits<double>::quiet_NaN());
    return output;
}

void check_output(const MetricOutput &output, int n_bits) {
    for (int bit = 0; bit < n_bits; ++bit) {
        test_assert(std::isfinite(output[bit]) && output[bit] >= 0 && output[bit] <= 1,
                    __FILE__, __LINE__, "Metric output was unwritten or outside [0, 1]");
    }
    test_assert(std::isnan(output[n_bits]), __FILE__, __LINE__, "Metric wrote past the requested bits");
}

void check_zero(const MetricOutput &output, int n_bits) {
    check_output(output, n_bits);
    for (int bit = 0; bit < n_bits; ++bit) {
        test_assert(output[bit] == 0, __FILE__, __LINE__, "Expected zero measured information");
    }
}

template <typename Function>
void rejects_empty_analysis(Function function) {
    bool rejected = false;
    try {
        function();
    } catch (const fixed_point_error &) {
        rejected = true;
    }
    test_assert(rejected, __FILE__, __LINE__, "Empty or null analysis input was accepted");
}

template <typename T>
void check_information_metrics(real_info_params<T> params) {
    const int n_bits = params.fixed_point ? 63
        : (params.mantissa_only ? std::numeric_limits<T>::digits - 1 : sizeof(T) * CHAR_BIT);
    const int varying_bit = params.fixed_point ? 62 : std::numeric_limits<T>::digits - 2;
    params.mask_zero = true;

    // Mixed pairs must use the same rows with masking enabled or disabled.
    T a[] = {T(0), T(1.5), T(1), T(1.5)};
    T b[] = {T(1.5), T(0), T(1.5), T(1)};
    auto masked_m = unwritten_output();
    auto masked_r = unwritten_output();
    mutual_information(params, a, b, 4, masked_m.data());
    redundancy(params, a, b, 4, masked_r.data());
    check_output(masked_m, n_bits);
    check_output(masked_r, n_bits);

    if (!params.fixed_point) {
        test_assert(double_equals(masked_m[varying_bit], 1) && double_equals(masked_r[varying_bit], 1),
                    __FILE__, __LINE__, "Complementary bits must have unit information and redundancy");
        auto h = unwritten_output();
        bit_count_entropy(params, a, 4, h.data());
        check_output(h, n_bits);
        test_assert(double_equals(h[varying_bit], 0.9182958340544896),
                    __FILE__, __LINE__, "Standalone entropy must still mask values independently");
    }

    params.mask_zero = false;
    auto unmasked_m = unwritten_output();
    auto unmasked_r = unwritten_output();
    mutual_information(params, a, b, 4, unmasked_m.data());
    redundancy(params, a, b, 4, unmasked_r.data());
    check_output(unmasked_m, n_bits);
    check_output(unmasked_r, n_bits);
    for (int bit = 0; bit < n_bits; ++bit) {
        test_assert(double_equals(masked_m[bit], unmasked_m[bit]) &&
                    double_equals(masked_r[bit], unmasked_r[bit]),
                    __FILE__, __LINE__, "Mixed pairs changed populations between paired metrics");
    }
    params.mask_zero = true;
    auto swapped_r = unwritten_output();
    redundancy(params, b, a, 4, swapped_r.data());
    check_output(swapped_r, n_bits);
    for (int bit = 0; bit < n_bits; ++bit) {
        test_assert(double_equals(swapped_r[bit], masked_r[bit]),
                    __FILE__, __LINE__, "Redundancy must be symmetric");
    }

    // A shifted view contains one masked 00 pair and four alternating 01/10 pairs.
    T shifted[] = {T(0), T(0), T(1.5), T(0), T(1.5), T(0)};
    std::array<size_t, 257> counts;
    counts.fill(19);
    size_t retained = 99;
    for (int repeat = 0; repeat < 2; ++repeat) {
        bit_pair_count_bits(params, shifted, shifted + 1, 5, counts.data(), &retained);
        test_assert(retained == 4, __FILE__, __LINE__, "Only the both-masked pair should be excluded");
        for (int bit = 0; bit < n_bits; ++bit) {
            const size_t *p = counts.data() + 4 * bit;
            test_assert(p[0] + p[1] + p[2] + p[3] == 4,
                        __FILE__, __LINE__, "Pair counts must overwrite previous output");
        }
        const size_t *p = counts.data() + 4 * varying_bit;
        test_assert(p[0] == 0 && p[1] == 2 && p[2] == 2 && p[3] == 0,
                    __FILE__, __LINE__, "Feature-boundary pair counts changed");
        test_assert(counts[4 * n_bits] == 19, __FILE__, __LINE__, "Pair counts exceeded their buffer");
    }
    auto spatial = unwritten_output();
    auto shifted_r = unwritten_output();
    bitwise_real_information(params, shifted, 6, spatial.data());
    redundancy(params, shifted, shifted + 1, 5, shifted_r.data());
    check_output(spatial, n_bits);
    check_output(shifted_r, n_bits);
    test_assert(double_equals(spatial[varying_bit], 1) && double_equals(shifted_r[varying_bit], 1),
                __FILE__, __LINE__, "Zero/nonzero transitions must contribute to every paired metric");

    // Disabling masking includes the extra 00 pair: probabilities are .2, .4, .4, 0.
    params.mask_zero = false;
    bitwise_real_information(params, shifted, 6, spatial.data());
    redundancy(params, shifted, shifted + 1, 5, shifted_r.data());
    const double marginal = -0.4 * std::log2(0.4) - 0.6 * std::log2(0.6);
    const double joint = -0.2 * std::log2(0.2) - 0.8 * std::log2(0.4);
    test_assert(double_equals(spatial[varying_bit], 2 * marginal - joint) &&
                double_equals(shifted_r[varying_bit], (2 * marginal - joint) / marginal),
                __FILE__, __LINE__, "Unmasked joint distribution gave the wrong metrics");
    params.mask_zero = true;

    // Masking uses a strict EPS threshold; a mixed pair retains the small member's bits.
    T threshold[] = {std::nextafter(T(EPS), T(0)), T(EPS), T(0)};
    T zeros[] = {T(0), T(0), T(0)};
    bit_pair_count_bits(params, threshold, zeros, 3, counts.data(), &retained);
    test_assert(retained == 1, __FILE__, __LINE__, "Masking must retain a value equal to EPS");
    T small_a[] = {T(EPS) / T(2), T(1.5)};
    T small_b[] = {T(1.5), T(EPS) / T(2)};
    redundancy(params, small_a, small_b, 2, masked_r.data());
    params.mask_zero = false;
    redundancy(params, small_a, small_b, 2, unmasked_r.data());
    for (int bit = 0; bit < n_bits; ++bit) {
        test_assert(double_equals(masked_r[bit], unmasked_r[bit]),
                    __FILE__, __LINE__, "A retained near-zero member must keep its actual bits");
    }
    params.mask_zero = true;

    // No retained samples: counts, entropy, mutual information, and redundancy are zero.
    T tiny[] = {T(0), T(EPS) / T(4), T(EPS) / T(2)};
    auto h = unwritten_output();
    auto m = unwritten_output();
    auto r = unwritten_output();
    bit_pair_count_bits(params, tiny, zeros, 3, counts.data(), &retained);
    test_assert(retained == 0, __FILE__, __LINE__, "All-masked pairs must leave no samples");
    for (int bit = 0; bit < 4 * n_bits; ++bit) {
        test_assert(counts[bit] == 0, __FILE__, __LINE__, "No-sample counts were not initialized");
    }
    bit_count_entropy(params, tiny, 3, h.data());
    mutual_information(params, tiny, zeros, 3, m.data());
    redundancy(params, tiny, zeros, 3, r.data());
    check_zero(h, n_bits);
    check_zero(m, n_bits);
    check_zero(r, n_bits);
    test_assert(preserved_information_template(params, tiny, tiny, 3) == 1 &&
                preserved_information_template(params, tiny, zeros, 3) == 0,
                __FILE__, __LINE__, "All-masked preservation must use bitwise identity");

    // Equal constant bits are lossless, unequal constants and a constant vs a varying bit are not.
    T constant[] = {T(1.5), T(1.5)};
    redundancy(params, zeros, constant, 2, r.data());
    check_output(r, n_bits);
    test_assert(r[varying_bit] == 0, __FILE__, __LINE__, "Different constant bits must have zero redundancy");
    redundancy(params, constant, constant, 2, r.data());
    for (int bit = 0; bit < n_bits; ++bit) {
        test_assert(r[bit] == 1, __FILE__, __LINE__, "Identical constant bits must have unit redundancy");
    }
    T varying[] = {T(0), T(1.5)};
    redundancy(params, varying, constant, 2, r.data());
    test_assert(r[varying_bit] == 0, __FILE__, __LINE__, "A constant bit cannot share a varying bit's information");

    T singleton = T(1.5);
    spatial = unwritten_output();
    bitwise_real_information(params, &singleton, 1, spatial.data());
    check_zero(spatial, n_bits);

    // Four adjacent pairs cover all outcomes equally; the reference is nonconstant but has zero MI.
    auto native = params;
    native.fixed_point = false;
    const int native_bits = native.mantissa_only ? std::numeric_limits<T>::digits - 1 : sizeof(T) * CHAR_BIT;
    T independent[] = {T(1), T(1.5), T(1.5), T(1), T(1)};
    T changed[] = {T(1.5), T(1.5), T(1.5), T(1), T(1)};
    spatial = unwritten_output();
    bitwise_real_information(native, independent, 5, spatial.data());
    check_zero(spatial, native_bits);
    test_assert(preserved_information_template(params, independent, independent, 5) == 1 &&
                preserved_information_template(params, independent, changed, 5) == 0,
                __FILE__, __LINE__, "Zero-information preservation must use bitwise identity");
    T signed_zeros[] = {T(0), -T(0), T(0)};
    test_assert(preserved_information_template(params, signed_zeros, signed_zeros, 3) == 1 &&
                preserved_information_template(params, signed_zeros, zeros, 3) == 0,
                __FILE__, __LINE__, "Constant preservation must distinguish signed-zero patterns");

    rejects_empty_analysis([&] { bit_count_entropy(params, a, 0, h.data()); });
    rejects_empty_analysis([&] { mutual_information(params, a, b, 0, m.data()); });
    rejects_empty_analysis([&] { redundancy(params, a, b, 0, r.data()); });
    rejects_empty_analysis([&] { bitwise_real_information(params, a, 0, spatial.data()); });
    rejects_empty_analysis([&] { bit_count_entropy<T>(params, nullptr, 1, h.data()); });
    rejects_empty_analysis([&] { mutual_information<T>(params, a, nullptr, 1, m.data()); });
    rejects_empty_analysis([&] { redundancy<T>(params, nullptr, b, 1, r.data()); });
    rejects_empty_analysis([&] { bitwise_real_information<T>(params, nullptr, 1, spatial.data()); });
}

template <typename T>
void check_modes() {
    std::string path = test_data_path("real_info.yaml");
    T domain[] = {T(0), T(1.5)};
    real_info_params<T> params(path.data());
    for (bool fixed : {false, true}) {
        params.fixed_point = fixed;
        if (fixed) params.fixed_range = get_fixed_point_init_params(domain, 2);
        for (bool mantissa : {false, true}) {
            params.mantissa_only = mantissa;
            check_information_metrics(params);
        }
    }
}

int main() {
    check_modes<float>();
    check_modes<double>();
    return 0;
}
