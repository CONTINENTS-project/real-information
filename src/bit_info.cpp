#include "bit_info.h"

namespace {
	thread_local std::string fixed_point_error_message;

	// Keep fixed-point validation exceptions inside the C++ boundary. Existing
	// C signatures stay unchanged; Python translates the message to ValueError.
	template <typename Function>
	auto fixed_point_call(Function function) -> std::invoke_result_t<Function> {
		using Result = std::invoke_result_t<Function>;
		fixed_point_error_message.clear();
		try {
			return function();
		} catch (const fixed_point_error &error) {
			fixed_point_error_message = error.what();
			if constexpr (std::is_void_v<Result>) {
				return;
			} else if constexpr (std::is_integral_v<Result>) {
				return -1;
			} else {
				return std::numeric_limits<Result>::quiet_NaN();
			}
		}
	}
}

extern "C" const char *real_info_last_fixed_point_error() {
	return fixed_point_error_message.empty() ? nullptr : fixed_point_error_message.c_str();
}

/**
 * Template instantiations for FORTRAN/C interface
 */

extern "C" {
	void shave_float(float *A, size_t n_elem, int n, float *shaved, char *config_file) {
		fixed_point_call([&] {
			if (n_elem == 0) return;
			real_info_params<float> params(config_file, A, n_elem);
			shave_bits<float>(params, A, n_elem, n, shaved);
		});
	}
}

extern "C" {
	void shave_double(double *A, size_t n_elem, int n, double *shaved, char *config_file) {
		fixed_point_call([&] {
			if (n_elem == 0) return;
			real_info_params<double> params(config_file, A, n_elem);
			shave_bits<double>(params, A, n_elem, n, shaved);
		});
	}
}

extern "C" {
	void halfshave_float(float *A, size_t n_elem, int n, float *halfshaved) {
		halfshave_template<float>(A, n_elem, n, halfshaved);
	}
}

extern "C" {
	void halfshave_double(double *A, size_t n_elem, int n, double *halfshaved) {
		halfshave_template<double>(A, n_elem, n, halfshaved);
	}
}

extern "C" {
	void bit_rounding_float(float *A, size_t n_elem, int n, float *bit_rounding_result) {
		bit_rounding_template<float>(A, n_elem, n, bit_rounding_result);
	}
}

extern "C" {
	void bit_rounding_double(double *A, size_t n_elem, int n, double *bit_rounding_result) {
		bit_rounding_template<double>(A, n_elem, n, bit_rounding_result);
	}
}

extern "C" {
	void set_float(float *A, size_t n_elem, int n, float *set) {
		set_template<float>(A, n_elem, n, set);
	}
}

extern "C" {
	void set_double(double *A, size_t n_elem, int n, double *set) {
		set_template<double>(A, n_elem, n, set);
	}
}

extern "C" {
	void groom_float(float *A, size_t n_elem, int n, float *groomed) {
		groom_template<float>(A, n_elem, n, groomed);
	}
}

extern "C" {
	void groom_double(double *A, size_t n_elem, int n, double *groomed) {
		groom_template<double>(A, n_elem, n, groomed);
	}
}

extern "C" {
	double preserved_information_float(float *A, float *B, size_t n_elem, char *config_file) {
		return fixed_point_call([&] {
			real_info_params<float> params(config_file, A, n_elem);
			return preserved_information_template<float>(params, A, B, n_elem);
		});
	}
}

extern "C" {
	double preserved_information_double(double *A, double *B, size_t n_elem, char *config_file) {
		return fixed_point_call([&] {
			real_info_params<double> params(config_file, A, n_elem);
			return preserved_information_template<double>(params, A, B, n_elem);
		});
	}
}

extern "C" {
	int pick_bits_to_shave_float(float *A, size_t n_elem, double tolerance, int nbits_old, char *config_file) {
		return fixed_point_call([&] {
			real_info_params<float> params(config_file, A, n_elem);
			return pick_bits_to_shave_template<float>(params, A, n_elem, tolerance, nbits_old);
		});
	}
}

extern "C" {
	int pick_bits_to_shave_double(double *A, size_t n_elem, double tolerance, int nbits_old, char *config_file) {
		return fixed_point_call([&] {
			real_info_params<double> params(config_file, A, n_elem);
			return pick_bits_to_shave_template<double>(params, A, n_elem, tolerance, nbits_old);
		});
	}
}

extern "C" {
	int pick_bits_to_shave_binary_search_float(float *A, size_t n_elem, double tolerance, int nbits_old, char *config_file , float *shaved) {
		return fixed_point_call([&] {
			real_info_params<float> params(config_file, A, n_elem);
			return pick_bits_to_shave_binary_search_template<float>(params, A, n_elem, tolerance, nbits_old, config_file, shaved);
		});
	}
}

extern "C" {
	int pick_bits_to_shave_binary_search_double(double *A, size_t n_elem, double tolerance, int nbits_old, char *config_file, double *shaved) {
		return fixed_point_call([&] {
			real_info_params<double> params(config_file, A, n_elem);
			return pick_bits_to_shave_binary_search_template<double>(params, A, n_elem, tolerance, nbits_old, config_file, shaved);
		});
	}
}

extern "C" {
	double entropy_double(double *A, size_t n_elem) {
		return bitpattern_entropy<double>(A, n_elem);
	}
}

extern "C" {
	double entropy_float(float *A, size_t n_elem) {
		return bitpattern_entropy<float>(A, n_elem);
	}
}
