#include "tests.h"
#include "bit_info.h"

int main() {
	float *data = new float[n_elements];
	const char *filename = "/work/d446/d446/ab_continents/cesm_workflow/real-information/tests/A.txt";
	read_input_data(filename, data, n_elements);

	size_t bit_counts[32];
	bit_count<float>(data, n_elements, bit_counts);

	size_t bit_counts_expected[32];
	read_yaml_list("/work/d446/d446/ab_continents/cesm_workflow/real-information/tests/expected.yaml", "bit_counts", bit_counts_expected, 32);
	for (int i = 0; i < 32; i++) {
		std::string message = "Bit count for bit " + std::to_string(i) + " does not match expected value " + std::to_string(bit_counts[i]) + " != " + std::to_string(bit_counts_expected[i]);
		test_assert(bit_counts[i] == bit_counts_expected[i], __FILE__, __LINE__, message);
	}

	bit_count_bits<float>(data, n_elements, bit_counts);
	for (int i = 0; i < 32; i++) {
		std::string message = "Bit count (bits class) for bit " + std::to_string(i) + " does not match expected value " + std::to_string(bit_counts[i]) + " != " + std::to_string(bit_counts_expected[i]);
		test_assert(bit_counts[i] == bit_counts_expected[i], __FILE__, __LINE__, message);
	}

	return 0;
}