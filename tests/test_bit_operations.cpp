#include "tests.h"
#include "bit_info.h"

#include <array>
#include <bit>
#include <limits>

template <typename Function>
void rejects_count(Function function) {
    bool rejected = false;
    try {
        function();
    } catch (const fixed_point_error &) {
        rejected = true;
    }
    test_assert(rejected, __FILE__, __LINE__, "Invalid bit count was accepted");
}

template <typename T>
using Word = std::conditional_t<std::is_same_v<T, float>, uint32_t, uint64_t>;

template <typename T>
using ArrayOperation = void (*)(const real_info_params<T> &, T *, size_t, int, T *);

template <typename T>
constexpr std::array<ArrayOperation<T>, 5> operations = {
    shave_bits<T>, set_bits<T>, halfshave_bits<T>, bit_rounding_bits<T>, groom_bits<T>
};

template <typename T>
constexpr std::array<void (*)(T *, size_t, int, T *), 5> legacy_operations = {
    shave_template<T>, set_template<T>, halfshave_template<T>, bit_rounding_template<T>, groom_template<T>
};

template <typename T>
constexpr std::array<void (bits<T>::*)(int), 4> member_operations = {
    &bits<T>::shave, &bits<T>::set, &bits<T>::halfshave, &bits<T>::round
};

template <typename T>
void check_native() {
    using UInt = Word<T>;
    constexpr int fraction = std::numeric_limits<T>::digits - 1;
    constexpr int width = sizeof(T) * CHAR_BIT;
    constexpr UInt sign = UInt{1} << (width - 1);
    constexpr UInt fraction_mask = (UInt{1} << fraction) - 1;
    constexpr UInt exponent = (std::numeric_limits<UInt>::max() >> 1) & ~fraction_mask;
    const UInt one = std::bit_cast<UInt>(T(1));
    const std::array<UInt, 12> patterns = {
        0, sign, 1, sign | 1, fraction_mask, fraction_mask + 1,
        one, one | sign, exponent, exponent | sign,
        exponent | (UInt{1} << (fraction - 1)) | 17, exponent | sign | 1
    };
    std::array<T, patterns.size()> input, output;
    for (size_t i = 0; i < input.size(); ++i) input[i] = std::bit_cast<T>(patterns[i]);
    std::string path = test_data_path("real_info.yaml");
    real_info_params<T> params(path.data());

    for (size_t op = 0; op < operations<T>.size(); ++op) {
        operations<T>[op](params, input.data(), input.size(), 0, output.data());
        test_assert(std::memcmp(input.data(), output.data(), sizeof(T) * input.size()) == 0,
                    __FILE__, __LINE__, "Zero-bit operation changed input bits");
        legacy_operations<T>[op](input.data(), input.size(), 0, output.data());
        test_assert(std::memcmp(input.data(), output.data(), sizeof(T) * input.size()) == 0,
                    __FILE__, __LINE__, "Legacy zero-bit operation changed input bits");
        auto in_place = input;
        operations<T>[op](params, in_place.data(), in_place.size(), 0, in_place.data());
        test_assert(std::memcmp(input.data(), in_place.data(), sizeof(T) * input.size()) == 0,
                    __FILE__, __LINE__, "In-place zero-bit operation changed input bits");
        operations<T>[op](params, nullptr, 0, fraction, nullptr);
        legacy_operations<T>[op](nullptr, 0, fraction, nullptr);

        for (int count : {-1, fraction + 1, width, std::numeric_limits<int>::min(), std::numeric_limits<int>::max()}) {
            output = input;
            rejects_count([&] { operations<T>[op](params, input.data(), input.size(), count, output.data()); });
            test_assert(std::memcmp(input.data(), output.data(), sizeof(T) * input.size()) == 0,
                        __FILE__, __LINE__, "Invalid operation changed output");
            rejects_count([&] { legacy_operations<T>[op](input.data(), input.size(), count, output.data()); });
        }
        for (int count : {1, fraction}) {
            operations<T>[op](params, input.data(), input.size(), count, output.data());
            std::array<T, patterns.size()> legacy_output;
            legacy_operations<T>[op](input.data(), input.size(), count, legacy_output.data());
            test_assert(std::memcmp(output.data(), legacy_output.data(), sizeof(T) * input.size()) == 0,
                        __FILE__, __LINE__, "Native and legacy operations disagree");
            for (size_t i = 8; i < patterns.size(); ++i) {
                test_assert(std::bit_cast<UInt>(output[i]) == patterns[i],
                            __FILE__, __LINE__, "Nonfinite payload changed");
            }
        }
    }
    for (auto operation : member_operations<T>) {
        bits<T> value(params, T(1));
        (value.*operation)(0);
        test_assert(std::bit_cast<UInt>(value.convert_back()) == one, __FILE__, __LINE__, "Member zero-bit identity failed");
        for (int count : {-1, fraction + 1, width}) {
            rejects_count([&] { (value.*operation)(count); });
            test_assert(std::bit_cast<UInt>(value.convert_back()) == one,
                        __FILE__, __LINE__, "Invalid member operation changed value");
        }
    }

    // Explicit IEEE patterns cover below/above midpoint and both tie parities.
    const std::array<std::pair<UInt, UInt>, 11> round_cases = {{
        {one + 7, one}, {one + 8, one}, {one + 9, one + 16},
        {one + 23, one + 16}, {one + 24, one + 32}, {one + 25, one + 32},
        {std::bit_cast<UInt>(T(2)) - 1, std::bit_cast<UInt>(T(2))},
        {8, 0}, {24, 32}, {fraction_mask, fraction_mask + 1}, {exponent - 1, exponent}
    }};
    for (UInt sign_bit : {UInt{0}, sign}) {
        for (const auto &[before, after] : round_cases) {
            T value = std::bit_cast<T>(before | sign_bit);
            T result;
            bit_rounding_bits(params, &value, 1, 4, &result);
            test_assert(std::bit_cast<UInt>(result) == (after | sign_bit),
                        __FILE__, __LINE__, "Nearest-even rounding pattern failed");
            bit_rounding_template(&value, 1, 4, &result);
            test_assert(std::bit_cast<UInt>(result) == (after | sign_bit),
                        __FILE__, __LINE__, "Legacy nearest-even rounding pattern failed");
            std::bitset<width> bitset(before | sign_bit);
            round_bitset(&bitset, 4);
            test_assert(bitset.to_ullong() == (after | sign_bit),
                        __FILE__, __LINE__, "Bitset nearest-even rounding pattern failed");
        }
        T value = std::bit_cast<T>(one | sign_bit);
        T result;
        halfshave_bits(params, &value, 1, 4, &result);
        test_assert(std::bit_cast<UInt>(result) == ((one + 8) | sign_bit),
                    __FILE__, __LINE__, "Half-shave did not set bit n-1");
        std::bitset<width> bitset(one | sign_bit);
        halfshave_bitset(&bitset, 4);
        test_assert(bitset.to_ullong() == ((one + 8) | sign_bit),
                    __FILE__, __LINE__, "Bitset half-shave used the wrong mask width");
        halfshave_bits(params, &value, 1, fraction, &result);
        test_assert(std::bit_cast<UInt>(result) == ((one + (UInt{1} << (fraction - 1))) | sign_bit),
                    __FILE__, __LINE__, "Maximum half-shave failed");
    }
    for (T value : {T(1.25), T(1.5), T(-1.25), T(-1.5)}) {
        T result;
        bit_rounding_bits(params, &value, 1, fraction, &result);
        const T expected = std::abs(value) < T(1.5) ? std::copysign(T(1), value) : std::copysign(T(2), value);
        test_assert(result == expected, __FILE__, __LINE__, "Maximum rounding count failed");
    }
    for (auto operation : {shave_bitset<width>, set_bitset<width>, halfshave_bitset<width>, round_bitset<width>}) {
        std::bitset<width> value(one + 3);
        operation(&value, 0);
        test_assert(value.to_ullong() == one + 3, __FILE__, __LINE__, "Bitset zero-bit identity failed");
        for (int count : {-1, fraction + 1, width}) rejects_count([&] { operation(&value, count); });
    }
}

template <typename T>
void check_fixed_point() {
    std::array<T, 5> input = {T(-1), T(-0.0), T(0), T(0.123456789), T(1)};
    std::array<T, 5> output;
    std::string path = test_data_path("fixed_point.yaml");
    real_info_params<T> params(path.data(), input.data(), input.size());
    for (auto operation : operations<T>) {
        operation(params, input.data(), input.size(), 0, output.data());
        test_assert(std::memcmp(input.data(), output.data(), input.size() * sizeof(T)) == 0,
                    __FILE__, __LINE__, "Fixed-point zero-bit identity lost bits");
        for (int count : {1, 31, 32, 52, 62, 63}) {
            operation(params, input.data(), input.size(), count, output.data());
            for (T value : output) {
                test_assert(std::isfinite(value) && value >= T(-1) && value <= T(1),
                            __FILE__, __LINE__, "Fixed-point output escaped its initialized range");
            }
        }
        for (int count : {-1, 64, std::numeric_limits<int>::min(), std::numeric_limits<int>::max()}) {
            output = input;
            rejects_count([&] { operation(params, input.data(), input.size(), count, output.data()); });
            test_assert(std::memcmp(input.data(), output.data(), input.size() * sizeof(T)) == 0,
                        __FILE__, __LINE__, "Invalid fixed-point operation changed output");
        }
    }
    for (const auto &[before, after] : std::array<std::pair<uint64_t, uint64_t>, 5>{{{7, 0}, {8, 0}, {9, 16}, {24, 32}, {25, 32}}}) {
        bits<T> value(params, T(-1));
        value.bit_set_64 = before;
        value.round(4);
        test_assert(value.bit_set_64 == after, __FILE__, __LINE__, "Fixed-point nearest-even tie failed");
    }
    for (int count : {1, 4, 31, 32, 52, 62, 63}) {
        bits<T> half(params, T(-1));
        half.halfshave(count);
        test_assert(half.bit_set_64 == (uint64_t{1} << (count - 1)),
                    __FILE__, __LINE__, "Fixed-point half-shave used native float width");
        for (auto operation : member_operations<T>) {
            bits<T> endpoint(params, T(1));
            (endpoint.*operation)(count);
            test_assert(endpoint.bit_set_64 <= static_cast<uint64_t>(FP<T>::max_code),
                        __FILE__, __LINE__, "Transformed code overflowed the upper endpoint");
            const T value = endpoint.convert_back();
            test_assert(std::isfinite(value) && value >= T(-1) && value <= T(1),
                        __FILE__, __LINE__, "Endpoint decoding failed");
        }
    }
    bits<T> upper(params, T(1));
    upper.round(63);
    test_assert(upper.convert_back() == T(1), __FILE__, __LINE__, "Rounded endpoint did not saturate");

    std::array<T, 3> constant = {T(2.5), T(2.5), T(2.5)};
    real_info_params<T> constant_params(path.data(), constant.data(), constant.size());
    for (auto operation : operations<T>) {
        operation(constant_params, constant.data(), constant.size(), 63, output.data());
        test_assert(std::equal(constant.begin(), constant.end(), output.begin()),
                    __FILE__, __LINE__, "Constant fixed-point data changed");
        rejects_count([&] { operation(constant_params, constant.data(), constant.size(), 64, output.data()); });
    }
}

void check_integer_rounding() {
    using real_info_detail::bit_operation;
    // Independent nearest-grid oracle: enumerate candidates and compare distances.
    for (uint32_t value = 0; value < 128; ++value) {
        for (int n = 1; n <= 5; ++n) {
            uint32_t best = 0;
            unsigned distance = std::numeric_limits<unsigned>::max();
            const uint32_t step = uint32_t{1} << n;
            for (uint32_t candidate = 0; candidate <= 160; candidate += step) {
                const unsigned delta = candidate > value ? candidate - value : value - candidate;
                if (delta < distance || (delta == distance && (candidate / step) % 2 == 0)) {
                    best = candidate;
                    distance = delta;
                }
            }
            test_assert(real_info_detail::transform_word(value, n, bit_operation::round) == best,
                        __FILE__, __LINE__, "32-bit integer rounding disagrees with nearest-grid oracle");
            test_assert(real_info_detail::transform_word(uint64_t(value), n, bit_operation::round) == best,
                        __FILE__, __LINE__, "64-bit integer rounding disagrees with nearest-grid oracle");
        }
    }
    rejects_count([] { real_info_detail::transform_word(uint64_t{1}, 64, bit_operation::round); });
    rejects_count([] { real_info_detail::transform_word(std::numeric_limits<uint64_t>::max(), 1, bit_operation::round); });
}

int main() {
    check_integer_rounding();
    check_native<float>();
    check_native<double>();
    check_fixed_point<float>();
    check_fixed_point<double>();
    return 0;
}
