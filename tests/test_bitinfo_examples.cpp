#include <iostream>
#include <fstream>
#include <random>
#include "../src/bit_info.h"

int main() {
	// Read data from file
	const size_t n_elem = 1000;
	float *data = new float[n_elem];
	std::ifstream infile("A.txt");
	for (size_t i = 0; i < n_elem; i++) {
		infile >> data[i];
	}

	// calculate bit pattern entropy
	double H = bitpattern_entropy<float>(data, n_elem);
	std::cout << "Bit pattern entropy: " << H << std::endl;

    // compute bit count
	size_t c[32];
	bit_count<float>(data, n_elem, c);
	std::cout << "Bit counts: ";
	for (int i = 0; i < 32; i++) {
		std::cout << c[i] << " ";
	}
	std::cout << "\n";

	bit_count_bits<float>(data, n_elem, c);
	std::cout << "Bit counts (bits class): ";
	for (int i = 0; i < 32; i++) {
		std::cout << c[i] << " ";
	}
	std::cout << "\n";

	// Compute bitcount entropy
	double *H_bits = new double[32];
	bit_count_entropy<float>(data, n_elem, H_bits);
	std::cout << "Bit count entropy: ";
	for (int i = 0; i < 32; i++) {
		std::cout << H_bits[i] << " ";
	}
	std::cout << "\n";

	// Compute mutual information between data and shaved version of data
	float *shaved_data = new float[n_elem];
	shave_template<float>(data, n_elem, 4, shaved_data);
	double *info = new double[32];
	mutual_information<float>(data, shaved_data, n_elem, info);
	std::cout << "Mutual information between data and shaved data: ";
	for (int i = 0; i < 32; i++) {
		std::cout << info[i] << " ";
	}
	std::cout << "\n";

	// Compute real information of sorted data
	float *sorted_data = new float[n_elem];
	std::memcpy(sorted_data, data, n_elem * sizeof(float));
	std::sort(sorted_data, sorted_data + n_elem);
	bitwise_real_information<float>(sorted_data, n_elem, info);
	std::cout << "Real information of sorted data: ";
	for (int i = 0; i < 32; i++) {
		std::cout << info[i] << " ";
	}
	std::cout << "\n";

	// Compute redundancy
	double *R = new double[32];
	redundancy<float>(data, shaved_data, n_elem, R);
	std::cout << "Redundancy between data and shaved data: ";
	for (int i = 0; i < 32; i++) {
		std::cout << R[i] << " ";
	}
	std::cout << "\n";

	// Compute preserved information
	double preserved_info = preserved_information_template<float>(data, shaved_data, n_elem);
	std::cout << "Preserved information in shaved data: " << preserved_info << "\n";

	return 0;
}