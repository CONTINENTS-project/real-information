#include "tests.h"
#include "bit_info.h"

int main() {
	float *data = new float[n_elements];
	const char *filename = "/work/d446/d446/ab_continents/cesm_workflow/real-information/tests/A.txt";
	read_input_data(filename, data, n_elements);

	double H[32];
	bit_count_entropy<float>(data, n_elements, H);

	double H_expected[32];
	read_yaml_list("/work/d446/d446/ab_continents/cesm_workflow/real-information/tests/expected.yaml", "bit_count_entropy", H_expected, 32);
	for (int i = 0; i < 32; i++) {
		std::string message = "Bit count entropy for bit " + std::to_string(i) + " does not match expected value " + std::to_string(H[i]) + " != " + std::to_string(H_expected[i]);
		test_assert(double_equals(H[i], H_expected[i]), __FILE__, __LINE__, message);
	}

	return 0;
}