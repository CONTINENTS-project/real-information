#include "tests.h"
#include "bit_info.h"

int main() {
	float *data = new float[n_elements];
	const std::string filename = test_data_path("A.txt");
	read_input_data(filename.c_str(), data, n_elements);

	float *shaved_data = new float[n_elements];
	shave_template<float>(data, n_elements, 4, shaved_data);
	std::string config_file = test_data_path("test_yaml_override.yaml");
	real_info_params<float> params(config_file.data(), data, n_elements);
	double preserved_information = preserved_information_template<float>(params, data, shaved_data, n_elements);

	double preserved_information_expected;
	read_yaml_value(test_data_path("expected.yaml"), "preserved_information", &preserved_information_expected);
	test_assert(double_equals(preserved_information, preserved_information_expected), __FILE__, __LINE__, "Preserved information does not match expected value");

	return 0;
}
