#include "tests.h"
#include "bit_info.h"

#include <filesystem>

extern "C" int pick_bits_to_shave_binary_search_double(
	double *, size_t, double, int, char *, double *);
extern "C" void shave_float(float *, size_t, int, float *, char *);

template <typename Function>
void expect_fixed_point_error(Function function, const char *message) {
	bool rejected = false;
	try {
		function();
	} catch (const fixed_point_error &) {
		rejected = true;
	}
	test_assert(rejected, __FILE__, __LINE__, message);
}

template <typename T>
void check_endpoints(T minimum, T maximum) {
	T input[] = {minimum, maximum};
	const auto range = get_fixed_point_init_params(input, 2);
	FP<T> lower(range, minimum);
	FP<T> upper(range, maximum);
	test_assert(lower.value == 0, __FILE__, __LINE__, "Minimum must map to zero");
	test_assert(upper.value == FP<T>::max_code, __FILE__, __LINE__, "Maximum must map to the safe upper code");
	test_assert(lower.convert_back() == minimum, __FILE__, __LINE__, "Minimum must decode exactly");
	test_assert(upper.convert_back() == maximum, __FILE__, __LINE__, "Maximum must decode exactly");
	test_assert(upper.convert_back_no_offset() == range.width, __FILE__, __LINE__, "Differences must remain in double precision");

	FP<T> encoded_zero(range);
	test_assert(encoded_zero.convert_back() == minimum, __FILE__, __LINE__, "Code construction must not encode floating-point zero");
	FP<T> invalid_code(range);
	invalid_code.value = std::numeric_limits<int64_t>::max();
	expect_fixed_point_error([&] { invalid_code.convert_back(); }, "Codes above the mapping must be rejected");
}

template <typename T>
void check_shaving_range(std::string &config, T minimum, T maximum) {
	T input[] = {minimum, maximum};
	T output[2];
	real_info_params<T> params(config.data(), input, 2);
	shave_bits(params, input, 2, 0, output);
	test_assert(std::memcmp(input, output, sizeof(input)) == 0,
		__FILE__, __LINE__, "Shaving zero bits must preserve the input exactly");
	shave_bits(params, input, 2, 63, output);
	test_assert(output[0] == minimum && output[1] == minimum,
		__FILE__, __LINE__, "Shaving all mapped bits must produce the lower endpoint");
}

template <typename T>
void check_constant(std::string &config, T constant) {
	T input[] = {constant, constant, constant};
	T output[3];
	real_info_params<T> params(config.data(), input, 3);
	test_assert(params.fixed_range.constant, __FILE__, __LINE__, "Constant mode must be explicit");
	test_assert(params.fixed_range.width == 0, __FILE__, __LINE__, "Constants must not require a scale");
	FP<T> encoded(params.fixed_range, constant);
	test_assert(encoded.value == 0 && encoded.convert_back() == constant,
		__FILE__, __LINE__, "Constant values must encode and decode without arithmetic");

	shave_bits(params, input, 3, 63, output);
	test_assert(std::memcmp(input, output, sizeof(input)) == 0,
		__FILE__, __LINE__, "Fixed-point shaving must preserve constants exactly");
	test_assert(pick_bits_to_shave_template(params, input, 3, 0.99, 0) == 0,
		__FILE__, __LINE__, "Constant selection must return zero bits shaved");
	test_assert(pick_bits_to_shave_binary_search_template(params, input, 3, 0.99, 0, config.data(), output) == 0,
		__FILE__, __LINE__, "Constant binary selection must return zero bits shaved");
	test_assert(std::memcmp(input, output, sizeof(input)) == 0,
		__FILE__, __LINE__, "Constant binary selection must fill its output");
	test_assert(preserved_information_template(params, input, output, 3) == 1.0,
		__FILE__, __LINE__, "Unchanged constants must be treated as lossless");
	output[0] = constant == T(0) ? T(1) : T(0);
	test_assert(preserved_information_template(params, input, output, 3) == 0.0,
		__FILE__, __LINE__, "Changed constants must not be treated as lossless");

	T singleton[] = {constant};
	real_info_params<T> singleton_params(config.data(), singleton, 1);
	test_assert(pick_bits_to_shave_template(singleton_params, singleton, 1, 0.99, 0) == 0,
		__FILE__, __LINE__, "Singleton conversion must use constant mode");
}

int main() {
	std::string config = (std::filesystem::path(__FILE__).parent_path() / "fixed_point.yaml").string();
	check_endpoints<float>(0.0f, 1.0f);
	check_endpoints<double>(0.0, 1.0);
	check_endpoints<float>(10.0f, 11.0f);
	check_endpoints<double>(10.0, 11.0);
	check_endpoints<float>(0.0f, 1.0e-30f);
	check_endpoints<float>(-std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
	check_shaving_range<float>(config, 0.0f, 1.0e-30f);
	check_shaving_range<float>(config, -std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
	check_shaving_range<double>(config, 10.0, 11.0);
	check_constant<float>(config, 0.0f);
	check_constant<double>(config, -2.5);

	float interior[] = {0.0f, 0.25f, 1.0f};
	FP<float> middle(get_fixed_point_init_params(interior, 3), interior[1]);
	test_assert(middle.convert_back() == interior[1], __FILE__, __LINE__, "Interior values must use normalized decoding");

	double zeros[] = {0.0, -0.0, 0.0};
	double output[3];
	real_info_params<double> zero_params(config.data(), zeros, 3);
	shave_bits(zero_params, zeros, 3, 63, output);
	test_assert(std::memcmp(zeros, output, sizeof(zeros)) == 0, __FILE__, __LINE__, "Constant shortcuts must preserve signed zeros");

	expect_fixed_point_error([&] { get_fixed_point_init_params(zeros, 0); }, "Empty initialization must be rejected");
	expect_fixed_point_error([] { get_fixed_point_init_params<double>(nullptr, 1); }, "Null input must be rejected");
	expect_fixed_point_error([&] { real_info_params<double> params(config.data(), zeros, 0); }, "Empty fixed-point parameters must be rejected");
	expect_fixed_point_error([&] { pick_bits_to_shave_template(zero_params, zeros, 0, 0.99, 0); }, "Empty selection must be rejected");
	expect_fixed_point_error([&] { preserved_information_template(zero_params, zeros, output, 0); }, "Empty preservation must be rejected");

	for (double invalid : {std::numeric_limits<double>::quiet_NaN(),
		std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()}) {
		double input[] = {0.0, invalid, 1.0};
		expect_fixed_point_error([&] { get_fixed_point_init_params(input, 3); }, "Nonfinite inputs must be rejected");
	}
	double too_wide[] = {-std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
	expect_fixed_point_error([&] { get_fixed_point_init_params(too_wide, 2); }, "Overflowing widths must be rejected");
	expect_fixed_point_error([&] { FP<double> invalid(get_fixed_point_init_params(zeros, 3), 1.0); }, "Values outside a constant range must be rejected");

	// Fixed-point errors must stay inside C++, with an inspectable C error.
	test_assert(pick_bits_to_shave_binary_search_double(zeros, 0, 0.99, 0, config.data(), output) == -1,
		__FILE__, __LINE__, "C selection must return an error for empty input");
	test_assert(real_info_last_fixed_point_error() != nullptr, __FILE__, __LINE__, "C validation errors must have a message");
	test_assert(pick_bits_to_shave_binary_search_double(zeros, 3, 0.99, 0, config.data(), output) == 0,
		__FILE__, __LINE__, "C selection must support constants");
	test_assert(real_info_last_fixed_point_error() == nullptr, __FILE__, __LINE__, "Successful calls must clear previous errors");
	test_assert(std::memcmp(zeros, output, sizeof(zeros)) == 0, __FILE__, __LINE__, "C selection must copy signed zeros");
	shave_float(nullptr, 0, 0, nullptr, nullptr);
	test_assert(real_info_last_fixed_point_error() == nullptr, __FILE__, __LINE__, "Empty elementwise shaving must be a no-op");
	return 0;
}
