#include <iostream>
#include <random>
#include <algorithm>
#include <array>
#include <bitset>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <omp.h>
#include "../src/bit_info.h"

int main() {
    size_t n_elements = 100;
    size_t n_elem = 1000000;
    size_t n_bits = 32;

    float *a = new float[n_elem];
    float *s = new float[n_elem];
    double *R = new double[n_bits];

    auto print_a = [&a]() {
        for (int i = 0; i < 10; i++) std::cout << a[i] << ' ';
        std::cout << '\n';
    }; 
    auto print_s = [&s]() {
        for (int i = 0; i < 10; i++) std::cout << s[i] << ' ';
        std::cout << '\n';
    }; 

    // Set random seed and generate random data
    int seed = 42; // You can choose any seed value
    std::mt19937 gen(seed); // Standard mersenne_twister_engine seeded with the seed
    //std::random_device rd;  // Will be used to obtain a seed for the random number engine
    //std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dist(0.0, 1.0);

    float period = n_elem;
    for (size_t i = 0; i < n_elem; i++) {
        a[i] = sin(i * 2 * M_PI / period) + 0.01 * dist(gen);
    }

    //float xs[n_elements];
    //double data[32];
    //shave_template(x, n_elements, 3, xs);
    //bit_count_entropy(x, n_elements, data);
    //redundancy(x, xs, n_elements, data);
    //bitwise_real_information(x, n_elements, data);
    //mutual_information(x, xs, n_elements, data);
    //auto p = preserved_information_template<float>(x, xs, n_elements);

    //float *a_shaved = new float[n_elem];
    //shave_template(a, n_elem, 3, a_shaved);
    //size_t *bit_pair_counts = new size_t[n_bits * 4];
    //bit_pair_count<float>(a, a_shaved, n_elem, bit_pair_counts);
    //for (int i = 0; i < n_bits; i++) {
    //    std::cout << "bit " << i << " pair count " << bit_pair_counts[i] << '\n';
    //}

    // Time the binary search for picking bits to shave
    auto start_time = std::chrono::high_resolution_clock::now();
    int shave_bits = pick_bits_to_shave_binary_search_template<float>(a, n_elem, 0.93, 0);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    int n_threads = omp_get_max_threads();
    std::cout << "Time taken for binary search to pick bits to shave: " << duration << " milliseconds\n";
    std::cout << "Bits to shave: " << shave_bits << '\n';
    std::cout << "Number of threads: " << n_threads << '\n';

    // Clean up dynamically allocated memory
    delete[] a;
    delete[] s;
    delete[] R;

    return 0;
}