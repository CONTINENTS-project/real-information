#include "tests.h"
#include "bit_info.h"

int main() {
	float *data = new float[n_elements];
	const std::string filename = test_data_path("A.txt");
	read_input_data(filename.c_str(), data, n_elements);

	double H[32];
	std::string config_file = test_data_path("real_info.yaml");
	real_info_params<float> params(config_file.data(), data, n_elements);
	bit_count_entropy<float>(params, data, n_elements, H);

	double H_expected[32];
	read_yaml_list(test_data_path("expected.yaml"), "bit_count_entropy", H_expected, 32);
	for (int i = 0; i < 32; i++) {
		std::string message = "Bit count entropy for bit " + std::to_string(i) + " does not match expected value " + std::to_string(H[i]) + " != " + std::to_string(H_expected[i]);
		test_assert(double_equals(H[i], H_expected[i]), __FILE__, __LINE__, message);
	}

	return 0;
}
