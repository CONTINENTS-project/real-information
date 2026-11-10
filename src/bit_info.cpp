#include "bit_info.h"

#include <cstdlib>

namespace {
	thread_local std::string fixed_point_error_message;

	// Report fixed-point validation failures and terminate the calling process.
	template <typename Function>
	auto fixed_point_call(Function function) -> std::invoke_result_t<Function> {
		fixed_point_error_message.clear();
		try {
			return function();
		} catch (const fixed_point_error &error) {
			fixed_point_error_message = error.what();
			std::cerr << "Fixed-point error: " << fixed_point_error_message << std::endl;
			std::exit(EXIT_FAILURE);
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
	void halfshave_float(float *A, size_t n_elem, int n, float *halfshaved, char *config_file) {
		fixed_point_call([&] {
			if (n_elem == 0) return;
			real_info_params<float> params(config_file, A, n_elem);
			halfshave_bits<float>(params, A, n_elem, n, halfshaved);
		});
	}
}

extern "C" {
	void halfshave_double(double *A, size_t n_elem, int n, double *halfshaved, char *config_file) {
		fixed_point_call([&] {
			if (n_elem == 0) return;
			real_info_params<double> params(config_file, A, n_elem);
			halfshave_bits<double>(params, A, n_elem, n, halfshaved);
		});
	}
}

extern "C" {
	void bit_rounding_float(float *A, size_t n_elem, int n, float *bit_rounding_result, char *config_file) {
		fixed_point_call([&] {
			if (n_elem == 0) return;
			real_info_params<float> params(config_file, A, n_elem);
			bit_rounding_bits<float>(params, A, n_elem, n, bit_rounding_result);
		});
	}
}

extern "C" {
	void bit_rounding_double(double *A, size_t n_elem, int n, double *bit_rounding_result, char *config_file) {
		fixed_point_call([&] {
			if (n_elem == 0) return;
			real_info_params<double> params(config_file, A, n_elem);
			bit_rounding_bits<double>(params, A, n_elem, n, bit_rounding_result);
		});
	}
}

extern "C" {
	void set_float(float *A, size_t n_elem, int n, float *set, char *config_file) {
		fixed_point_call([&] {
			if (n_elem == 0) return;
			real_info_params<float> params(config_file, A, n_elem);
			set_bits<float>(params, A, n_elem, n, set);
		});
	}
}

extern "C" {
	void set_double(double *A, size_t n_elem, int n, double *set, char *config_file) {
		fixed_point_call([&] {
			if (n_elem == 0) return;
			real_info_params<double> params(config_file, A, n_elem);
			set_bits<double>(params, A, n_elem, n, set);
		});
	}
}

extern "C" {
	void groom_float(float *A, size_t n_elem, int n, float *groomed, char *config_file) {
		fixed_point_call([&] {
			if (n_elem == 0) return;
			real_info_params<float> params(config_file, A, n_elem);
			groom_bits<float>(params, A, n_elem, n, groomed);
		});
	}
}

extern "C" {
	void groom_double(double *A, size_t n_elem, int n, double *groomed, char *config_file) {
		fixed_point_call([&] {
			if (n_elem == 0) return;
			real_info_params<double> params(config_file, A, n_elem);
			groom_bits<double>(params, A, n_elem, n, groomed);
		});
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

// New symbols keep the legacy C/Fortran ABI intact. These entry points report
// errors to ctypes instead of terminating the Python worker on invalid caps.
template <typename T>
int capped_selection(T *input, size_t count, double tolerance, int previous,
                     char *config, int cap, T *output) {
    fixed_point_error_message.clear();
    try {
        if (cap < 0) throw fixed_point_error("max_shave_bits must be nonnegative.");
        real_info_params<T> params(config, input, count);
        return pick_bits_to_shave_binary_search_capped_template(
            params, input, count, tolerance, previous, cap, output);
    } catch (const std::exception &error) {
        fixed_point_error_message = error.what();
        return -1;
    }
}

extern "C" int pick_bits_to_shave_binary_search_capped_float(
    float *input, size_t count, double tolerance, int previous, char *config, int cap, float *output) {
    return capped_selection(input, count, tolerance, previous, config, cap, output);
}

extern "C" int pick_bits_to_shave_binary_search_capped_double(
    double *input, size_t count, double tolerance, int previous, char *config, int cap, double *output) {
    return capped_selection(input, count, tolerance, previous, config, cap, output);
}

extern "C" {
	double entropy_float(float *A, size_t n_elem) {
		return bitpattern_entropy<float>(A, n_elem);
	}
}
