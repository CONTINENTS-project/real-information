#include "tests.h"
#include "bit_info.h"

int main() {
	float *data = new float[n_elements];
	const char *filename = "/work/d446/d446/ab_continents/cesm_workflow/real-information/tests/A.txt";
	read_input_data(filename, data, n_elements);

	float *shaved_data = new float[n_elements];
	shave_template<float>(data, n_elements, 4, shaved_data);
	double *info = new double[32];
	mutual_information<float>(data, shaved_data, n_elements, info);

	double info_expected[32];
	read_yaml_list("/work/d446/d446/ab_continents/cesm_workflow/real-information/tests/expected.yaml", "mutual_information", info_expected, 32);
	for (int i = 0; i < 32; i++) {
		std::string message = "Mutual information for bit " + std::to_string(i) + " does not match expected value " + std::to_string(info[i]) + " != " + std::to_string(info_expected[i]);
		test_assert(double_equals(info[i], info_expected[i]), __FILE__, __LINE__, message);
	}

	return 0;
}