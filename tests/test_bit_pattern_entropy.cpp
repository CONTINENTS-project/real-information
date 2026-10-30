#include "tests.h"
#include "bit_info.h"

int main() {
	float *data = new float[n_elements];
	const std::string filename = test_data_path("A.txt");
	read_input_data(filename.c_str(), data, n_elements);

	double H = bitpattern_entropy<float>(data, n_elements);
	double H_expected;
	read_yaml_value(test_data_path("expected.yaml"), "bit_pattern_entropy", &H_expected);
	std::string message = "Bit pattern entropy does not match expected value " + std::to_string(H) + " != " + std::to_string(H_expected);
	test_assert(float_equals(H, H_expected), __FILE__, __LINE__, message);

	return 0;
}
