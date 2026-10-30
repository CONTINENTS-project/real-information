#include "tests.h"
#include "bit_info.h"

int main() {
	float *data = new float[n_elements];
	const std::string filename = test_data_path("A.txt");
	read_input_data(filename.c_str(), data, n_elements);

	std::string config_file = test_data_path("real_info.yaml");
	real_info_params<float> params(config_file.data(), data, n_elements);

	size_t bit_counts[32];
	size_t n_non_zero = 0;
	bit_count<float>(params, data, n_elements, bit_counts, &n_non_zero);

	size_t bit_counts_expected[32];
	read_yaml_list(test_data_path("expected.yaml"), "bit_counts", bit_counts_expected, 32);
	for (int i = 0; i < 32; i++) {
		std::string message = "Bit count for bit " + std::to_string(i) + " does not match expected value " + std::to_string(bit_counts[i]) + " != " + std::to_string(bit_counts_expected[i]);
		test_assert(bit_counts[i] == bit_counts_expected[i], __FILE__, __LINE__, message);
	}

	bit_count_bits<float>(params, data, n_elements, bit_counts, &n_non_zero);
	for (int i = 0; i < 32; i++) {
		std::string message = "Bit count (bits class) for bit " + std::to_string(i) + " does not match expected value " + std::to_string(bit_counts[i]) + " != " + std::to_string(bit_counts_expected[i]);
		test_assert(bit_counts[i] == bit_counts_expected[i], __FILE__, __LINE__, message);
	}

	return 0;
}
