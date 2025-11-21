#include "bit_info.h"

/**
 * Template instantiations for FORTRAN/C interface
 */

extern "C" {
	void shave_float(float *A, size_t n_elem, int n, float *shaved) {
		shave_template<float>(A, n_elem, n, shaved);
	}
}

extern "C" {
	void shave_double(double *A, size_t n_elem, int n, double *shaved) {
		shave_template<double>(A, n_elem, n, shaved);
	}
}

extern "C" {
	double preserved_information_float(float *A, float *B, size_t n_elem) {
		return preserved_information_template<float>(A, B, n_elem);
	}
}

extern "C" {
	double preserved_information_double(double *A, double *B, size_t n_elem) {
		return preserved_information_template<double>(A, B, n_elem);
	}
}

extern "C" {
	int pick_bits_to_shave_float(float *A, size_t n_elem, double tolerance, int nbits_old) {
		return pick_bits_to_shave_template<float>(A, n_elem, tolerance, nbits_old);
	}
}

extern "C" {
	int pick_bits_to_shave_double(double *A, size_t n_elem, double tolerance, int nbits_old) {
		return pick_bits_to_shave_template<double>(A, n_elem, tolerance, nbits_old);
	}
}

extern "C" {
	int pick_bits_to_shave_binary_search_float(float *A, size_t n_elem, double tolerance, int nbits_old) {
		return pick_bits_to_shave_binary_search_template<float>(A, n_elem, tolerance, nbits_old);
	}
}

extern "C" {
	int pick_bits_to_shave_binary_search_double(double *A, size_t n_elem, double tolerance, int nbits_old) {
		return pick_bits_to_shave_binary_search_template<double>(A, n_elem, tolerance, nbits_old);
	}
}