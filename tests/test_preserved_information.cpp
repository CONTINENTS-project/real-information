#include "tests.h"
#include "bit_info.h"

int main() {
	float *data = new float[n_elements];
	const char *filename = "/work/d446/d446/ab_continents/cesm_workflow/real-information/tests/A.txt";
	read_input_data(filename, data, n_elements);

	float *shaved_data = new float[n_elements];
	shave_template<float>(data, n_elements, 4, shaved_data);
	double preserved_information = preserved_information_template<float>(data, shaved_data, n_elements);

	double preserved_information_expected;
	read_yaml_value("/work/d446/d446/ab_continents/cesm_workflow/real-information/tests/expected.yaml", "preserved_information", &preserved_information_expected);
	test_assert(double_equals(preserved_information, preserved_information_expected), __FILE__, __LINE__, "Preserved information does not match expected value");

	return 0;
}