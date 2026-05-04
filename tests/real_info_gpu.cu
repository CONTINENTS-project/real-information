#include <iostream>
#include <cstdio>
#include "../src/cuda_bitset.h"

__global__ void shift_left(float *a, float *c, bits<float> *bits_a, bits<float> *bits_b, size_t n_elem) {
	size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < n_elem) {
		bits_a[idx] = bits<float>(a[idx]);
		auto tmp = bits_a[idx];
		tmp <<= 1;
		bits_b[idx] = tmp;

		c[idx] = bits_b[idx].convert_back();
		printf("device i %lu %f\n", idx, c[idx]);
	}
}

int main() {
	size_t n_elem = 10;
	float *a_h, *b_h, *c_h;
	float *a_d, *b_d, *c_d;
	bits<float> *bitset_a_d, *bitset_b_d;

    a_h = new float[n_elem];
	b_h = new float[n_elem];
	c_h = new float[n_elem];
	cudaMalloc(&bitset_a_d, sizeof(bits<float>) * n_elem);
	cudaMalloc(&bitset_b_d, sizeof(bits<float>) * n_elem);
	cudaMalloc(&a_d, sizeof(float) * n_elem);
	cudaMalloc(&b_d, sizeof(float) * n_elem);
	cudaMalloc(&c_d, sizeof(float) * n_elem);

	for (size_t i = 0; i < n_elem; i++) {
		a_h[i] = i;
	}
	for (size_t i = 0; i < n_elem; i++) {
		std::cout << "start host i " << i << " " << a_h[i] << std::endl;
	}

	cudaMemcpy(a_d, a_h, sizeof(float) * n_elem, cudaMemcpyHostToDevice);

	shift_left<<<(n_elem + 255) / 256, 256>>>(a_d, c_d, bitset_a_d, bitset_b_d, n_elem);

	cudaMemcpy(c_h, c_d, sizeof(float) * n_elem, cudaMemcpyDeviceToHost);

	for (size_t i = 0; i < n_elem; i++) {
		std::cout << "end host i " << i << " " << c_h[i] << std::endl;
	}

	return 0;
}