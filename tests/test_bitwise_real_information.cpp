#include "tests.h"
#include "bit_info.h"

int main() {
	float *data = new float[n_elements];
	const std::string filename = test_data_path("A.txt");
	read_input_data(filename.c_str(), data, n_elements);

	float *sorted_data = new float[n_elements];
	std::memcpy(sorted_data, data, n_elements * sizeof(float));
	std::sort(sorted_data, sorted_data + n_elements);
	double *info = new double[32];
	std::string config_file = test_data_path("real_info.yaml");
	real_info_params<float> params(config_file.data(), sorted_data, n_elements);
	bitwise_real_information<float>(params, sorted_data, n_elements, info);

	double info_expected[32];
	read_yaml_list(test_data_path("expected.yaml"), "bitwise_real_information", info_expected, 32);
	for (int i = 0; i < 32; i++) {
		std::string message = "Bitwise real information for bit " + std::to_string(i) + " does not match expected value " + std::to_string(info[i]) + " != " + std::to_string(info_expected[i]);
		test_assert(double_equals(info[i], info_expected[i]), __FILE__, __LINE__, message);
	}

	return 0;
}
