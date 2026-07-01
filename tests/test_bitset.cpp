#include <iostream>
#include <random>
#include <chrono>
#include <bitset>
#include "../src/bit_info.h"

#define N_ELEMENTS 1000

int main() {
	//int p[4];
	//float a_f = 1.61e-43;
	//float b_f = 1.19e-43;
	//bits<float> a(1.61e-43);
	//bits<float> b(1.19e-43);
	//auto a_bitset = reinterpret_cast<typename HELP<float>::type *>(&a_f);
	//auto b_bitset = reinterpret_cast<typename HELP<float>::type *>(&b_f);

	//for (int j = 0; j < 7; j++) {
	//	bit_pair_count_bits(a, b, j, p);
	//	//std::cout << "bit position        " << j << ": 00=" << p[0] << " 01=" << p[1] << " 10=" << p[2] << " 11=" << p[3] << "\n";
	//	bit_pair_count<32>(*a_bitset, *b_bitset, j, p);
	//	//std::cout << "bitset bit position " << j << ": 00=" << p[0] << " 01=" << p[1] << " 10=" << p[2] << " 11=" << p[3] << "\n";
	//}

	float A[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
	float B[10] = {10, 9, 8, 7, 6, 5, 4, 3, 2, 1};

    const int n_joint_counts = 4;
	size_t pair_counts[32 * n_joint_counts];
	size_t pair_counts_bits[32 * n_joint_counts];

    for (int i = 0; i < 32 * n_joint_counts; i++) {
		pair_counts[i] = 0;
		pair_counts_bits[i] = 0;
	}

    bit_pair_count_bits(A, B, 10, pair_counts_bits);
	bit_pair_count(A, B, 10, pair_counts);

	for (int i = 0; i < 32 * n_joint_counts; i++) {
		if (pair_counts[i] != pair_counts_bits[i]) {
			std::cout << "Mismatch at index " << i << ": " << pair_counts[i] << " vs " << pair_counts_bits[i] << "\n";
		}
	}
	std::cout << "Test completed.\n";

	size_t bit_counts[32];
	size_t bit_counts_bits[32];

	bit_count_bits(A, 10, bit_counts_bits);
	bit_count(A, 10, bit_counts);
	for (int i = 0; i < 32; i++) {
		if (bit_counts[i] != bit_counts_bits[i]) {
			std::cout << "Mismatch at index " << i << ": " << bit_counts[i] << " vs " << bit_counts_bits[i] << "\n";
		}
	}
	std::cout << "Bit count test completed.\n";

	// Create a random array of floats
	float random_floats[N_ELEMENTS];
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0.0, 1.0);
	for (int i = 0; i < N_ELEMENTS; i++) {
		random_floats[i] = dis(gen);
	}

	float shaved_floats_bits[N_ELEMENTS];
	float shaved_floats_bitset[N_ELEMENTS];

	for (int n = 0; n < 31; n++) {
		set_bits(random_floats, N_ELEMENTS, n, shaved_floats_bits);
		set_template(random_floats, N_ELEMENTS, n, shaved_floats_bitset);

		for (int i = 0; i < N_ELEMENTS; i++) {
			if (shaved_floats_bits[i] != shaved_floats_bitset[i]) {
				std::cout << "Mismatch at index " << i << " for n=" << n << ": " << shaved_floats_bits[i] << " vs " << shaved_floats_bitset[i] << "\n";
			}
		}
	}
	std::cout << "Shave test completed.\n";
	
	// Create a random array of doubles
	double random_doubles[N_ELEMENTS];
	for (int i = 0; i < N_ELEMENTS; i++) {
		random_doubles[i] = dis(gen);
	}

	double shaved_doubles_bits[N_ELEMENTS];
	double shaved_doubles_bitset[N_ELEMENTS];

	for (int n = 0; n < 31; n++) {
		set_bits(random_doubles, N_ELEMENTS, n, shaved_doubles_bits);
		set_template(random_doubles, N_ELEMENTS, n, shaved_doubles_bitset);

		for (int i = 0; i < N_ELEMENTS; i++) {
			if (shaved_doubles_bits[i] != shaved_doubles_bitset[i]) {
				std::cout << "Mismatch at index " << i << " for n=" << n << ": " << shaved_doubles_bits[i] << " vs " << shaved_doubles_bitset[i] << "\n";
			}
		}
	}
	std::cout << "Shave test completed.\n";

	// add timer
	auto start = std::chrono::high_resolution_clock::now();
	int n_bits = pick_bits_to_shave_binary_search_template(random_floats, N_ELEMENTS, 0.99, 0);
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
	std::cout << "Optimal number of bits to shave: " << n_bits << "\n";
	std::cout << "Time taken: " << duration.count() << " microseconds\n";

	return 0;
}