#ifndef BIT_INFO_H
#define BIT_INFO_H

#include <algorithm>
#include <concepts>
#include <climits>
#include <bitset>
#include <iostream>
#include <cmath>
#include <cstring>
#include <mpi.h>

inline int get_n_exponent_bits(int nbits);
const bool mantissa_only = true;

void get_start_end_bits(int n_bits, int *start, int *end) {
	if (mantissa_only) {
		if (n_bits == 32) {
			*end = 23;
		} else if (n_bits == 64) {
			*end = 52;
		}
	} else {
		*end = n_bits;
	}
	*start = 0;
}

template<typename T>
/**
 * @brief Utility wrapper for bit-level operations on floating-point values.
 *
 * The class stores the raw IEEE-754 bit pattern of a float or double and
 * exposes helpers to shave, set, halfshave, and round least-significant bits.
 *
 * @tparam T Floating-point type (`float` or `double`).
 */
class bits{
	public:
		uint32_t bit_set_32;
		uint64_t bit_set_64;

		/**
		 * @brief Construct from a floating-point value by copying its bit pattern.
		 * @param x Input value.
		 */
		bits(T x) {
			if (std::is_same_v<T, float>) {
				std::memcpy(&bit_set_32, &x, sizeof(float));
			} else if (std::is_same_v<T, double>) {
				std::memcpy(&bit_set_64, &x, sizeof(double));
			}
		}

		/**
		 * @brief Reconstruct the floating-point value from the stored bits.
		 * @return Value of type `T` represented by the current bit pattern.
		 */
		T convert_back() {
			T result;
			if (std::is_same_v<T, float>) {
				std::memcpy(&result, &bit_set_32, sizeof(float));
			} else if (std::is_same_v<T, double>) {
				std::memcpy(&result, &bit_set_64, sizeof(double));
			}
			return result;
		}

		/**
		 * @brief Shave the lowest `n` bits.
		 * @param n Number of least-significant bits to shave.
		 */
		void shave(int n) {
			if (std::is_same_v<T, float>) {
				uint32_t mask = 0xFFFFFFFF << n;
				bit_set_32 &= mask;
			} else if (std::is_same_v<T, double>) {
				uint64_t mask = 0xFFFFFFFFFFFFFFFF << n;
				bit_set_64 &= mask;
			}
		}

		/**
		 * @brief Set the lowest `n` bits to 1.
		 * @param n Number of least-significant bits to set.
		 */
		void set(int n) {
			if (std::is_same_v<T, float>) {
				uint32_t mask = 0xFFFFFFFF << n;
				bit_set_32 |= ~mask;
			} else if (std::is_same_v<T, double>) {
				uint64_t mask = 0xFFFFFFFFFFFFFFFF << n;
				bit_set_64 |= ~mask;
			}
		}

		/**
		 * @brief Shave the lowest `n` bits and set the highest cleared bit.
		 * @param n Number of least-significant bits to modify.
		 */
		void halfshave(int n) {
			int nbits = sizeof(T) * CHAR_BIT;
			if (std::is_same_v<T, float>) {
				uint32_t keep_mask = 0xFFFFFFFF << n;
				bit_set_32 &= keep_mask;
				uint32_t round_up_mask = 0x80000000 >> (nbits - n - 1);
				bit_set_32 |= round_up_mask;
			} else if (std::is_same_v<T, double>) {
				uint64_t keep_mask = 0xFFFFFFFFFFFFFFFF << n;
				bit_set_64 &= keep_mask;
				uint64_t round_up_mask = 0x8000000000000000 >> (nbits - n - 1);
				bit_set_64 |= round_up_mask;
			}
		}

		/**
		 * @brief Bitwise round the `n` least-significant bits
		 * @param n Number of least-significant bits to round away.
		 */
		void round(int n) {
			int nbits = sizeof(T) * CHAR_BIT;
			int mantissa_bits = nbits - get_n_exponent_bits(nbits) - 1;
			int keepbits = mantissa_bits - n;
			int shift = mantissa_bits - keepbits - 1;
			shift -= (keepbits == mantissa_bits) ? 1 : 0;

			if (std::is_same_v<T, float>) {
				uint32_t upper_mask;
				if (nbits - n + 1 >= 32) {
					upper_mask = 0x00000000; // bitshift leads to undefined behaviour if greater than number of bits, and what we really want are all zeros
				} else {
				 	upper_mask = 0xffffffff >> (nbits - n + 1);
				}
				uint32_t lower_mask = 0xFFFFFFFF << n;
				uint32_t temp_bit_set_32 = bit_set_32 + upper_mask + ((bit_set_32 >> shift) & 1);
				bit_set_32 = temp_bit_set_32 & lower_mask;
			} else if (std::is_same_v<T, double>) {
				uint64_t upper_mask;
				if (nbits - n + 1 >= 64) {
					upper_mask = 0x0000000000000000; // bitshift leads to undefined behaviour if greater than number of bits, and what we really want are all zeros
				} else {
				 	upper_mask = 0xFFFFFFFFFFFFFFFF >> (nbits - n + 1);
				}
				uint64_t lower_mask = 0xFFFFFFFFFFFFFFFF << n;
				uint64_t temp_bit_set_64 = bit_set_64 + upper_mask + ((bit_set_64 >> shift) & 1);
				bit_set_64 = temp_bit_set_64 & lower_mask;
			}
		}

};

/** 
 * @brief Shave the lowest `n` bits from each element in an array.
 * @param A Input array of type `T`.
 * @param n_elem Number of elements in the array.
 * @param n Number of least-significant bits to shave.
 * @param shaved Output array of shaved values.
 */
template <typename T>
void shave_bits(T *A, size_t n_elem, int n, T *shaved) {
	for (size_t i = 0; i < n_elem; i++) {
		bits<T> b(A[i]);
		b.shave(n);
		shaved[i] = b.convert_back();
	}
}

/** 
 * @brief Half-shave the lowest `n` bits from each element in an array.
 * @param A Input array of type `T`.
 * @param n_elem Number of elements in the array.
 * @param n Number of least-significant bits to half-shave.
 * @param halfshaved Output array of half-shaved values.
 */
template <typename T>
void halfshave_bits(T *A, size_t n_elem, int n, T *halfshaved) {
	for (size_t i = 0; i < n_elem; i++) {
		bits<T> b(A[i]);
		b.halfshave(n);
		halfshaved[i] = b.convert_back();
	}
}

/** 
 * @brief Bitwise round the `n` least-significant bits from each element in an array.
 * @param A Input array of type `T`.
 * @param n_elem Number of elements in the array.
 * @param n Number of least-significant bits to round.
 * @param bit_rounding_result Output array of rounded values.
 */
template <typename T>
void bit_rounding_bits(T *A, size_t n_elem, int n, T *bit_rounding_result) {
	for (size_t i = 0; i < n_elem; i++) {
		bits<T> b(A[i]);
		b.round(n);
		bit_rounding_result[i] = b.convert_back();
	}
}

/** 
 * @brief Set the lowest `n` bits for each element in an array.
 * @param A Input array of type `T`.
 * @param n_elem Number of elements in the array.
 * @param n Number of least-significant bits to set.
 * @param set Output array of set values.
 */
template <typename T>
void set_bits(T *A, size_t n_elem, int n, T *set) {
	for (size_t i = 0; i < n_elem; i++) {
		bits<T> b(A[i]);
		b.set(n);
		set[i] = b.convert_back();
	}
}

static int bit_pair_count_calls = 0;

template<std::floating_point T>
struct HELP {
	using type = std::bitset<sizeof(T) * CHAR_BIT>;
};

static int n_exponent_bits_float = 8;
static int n_exponent_bits_double = 11;
static int n_joint_counts = 4;
static int n_counts = 2;

/** 
 *  @brief Gets the number of exponent bits for a given bit size following IEEE 754 standard.
 *  @param nbits The bit size.
 *  @return The number of exponent bits.
 */
inline int get_n_exponent_bits(int nbits) {
	if (nbits == 32) {
		return 23;
	} else if (nbits == 64) {
		return 52;
	} else {
		std::cerr << "Unsupported bit size: " << nbits << "\n";
		return 0;
	}
}

// XXX: temporarily dump these helper functions here
/**
 * @brief Counts the number of bits set to 1 in a specific bit position across an array.
 * @param A Input array of type `T`.
 * @param n_elem Number of elements in the array.
 * @param i Bit position to count (0 for least significant bit).
 */
template <typename T>
size_t bit_count(T *A, size_t n_elem, int i) {
	size_t count = 0;

	for (size_t j = 0; j < n_elem; j++) {
		auto *b0 = reinterpret_cast<typename HELP<T>::type*>(&A[j]);
		count += (*b0)[i];
	}

	return count;
}

/**
 * @brief Counts the number of bits set to 1 in a specific bit position across an array using the `bits` class.
 * @param A Input array of type `T`.
 * @param n_elem Number of elements in the array.
 * @param i Bit position to count (0 for least significant bit).
 * @return The count of bits set to 1 in the specified bit position across the array.
 */
template <typename T>
size_t bit_count_bits(T *A, size_t n_elem, int i) {
	size_t count = 0;
	for (size_t j = 0; j < n_elem; j++) {
		bits<T> b(A[j]);
		if (std::is_same_v<T, float>) {
			uint32_t mask = 0x00000001;
			mask <<= i;
			count += (b.bit_set_32 & mask) ? 1 : 0;
		} else if (std::is_same_v<T, double>) {
			uint64_t mask = 0x0000000000000001;
			mask <<= i;
			count += (b.bit_set_64 & mask) ? 1 : 0;
		}
	}

	return count;
}

/**
 * @brief Counts the number of bits set to 1 in each bit position across an array.
 * @param A Input array of type `T`.
 * @param n_elem Number of elements in the array.
 * @param c Output array to store the count of bits set to 1 for each bit position.
 */
template <typename T>
void bit_count(T *A, size_t n_elem, size_t *c) {
	int n_bits = sizeof(A[0]) * CHAR_BIT;
	int start_bit, end_bit;
	get_start_end_bits(n_bits, &start_bit, &end_bit);
	for (int i = start_bit; i < end_bit; i++) {
		c[i - start_bit] = bit_count(A, n_elem, i);
	}
}
template <typename T>
void bit_count_bits(T *A, size_t n_elem, size_t *c) {
	int n_bits = sizeof(A[0]) * CHAR_BIT;
	int start_bit, end_bit;
	get_start_end_bits(n_bits, &start_bit, &end_bit);
	for (int i = start_bit; i < end_bit; i++) {
		c[i - start_bit] = bit_count_bits(A, n_elem, i);
	}
}

/**
 * @brief Computes the Shannon entropy of a probability distribution.
 * @param A Input array of probabilities (must sum to 1).
 * @param n_elem Number of elements in the array.
 * @return The Shannon entropy of the distribution represented by A.
 */
double entropy(double *A, size_t n_elem) {
	double H = 0;

	for (size_t i = 0; i < n_elem; i++) {
		if (A[i] > 0) {
			H += -A[i] * std::log2(A[i]);
		}
	}

	return H;
}

/** 
 * @brief Computes the entropy of a bit pattern distribution.
 * @param A Input array of type `T`.
 * @param n_elem Number of elements in the array.
 * @return The entropy of the bit pattern distribution.
 */
template <typename T>
double bitpattern_entropy(T *A, size_t n_elem) {
  double p;
  double H = 0;
  size_t count = 1;

  T *A_sorted = new T[n_elem];
  std::memcpy(A_sorted, A, n_elem * sizeof(T));
  std::sort(A_sorted, A_sorted + n_elem);
  for (size_t i = 0; i < n_elem-1; i++) {
    if (A_sorted[i+1] == A_sorted[i]) {
      count += 1;
    } else {
      p = static_cast<double>(count)/static_cast<double>(n_elem);
      H += -p * std::log2(p);
      count = 1;
    }
  }

  p = static_cast<double>(count)/static_cast<double>(n_elem);
  H += -p * std::log2(p);
  delete[] A_sorted;
  return H;
}

/**
 * @brief Helper function to compute pair counts at a particular position in two numbers represented as bits.
 * @param a First number represented as bits.
 * @param b Second number represented as bits.
 * @param j Bit position to compare.
 * @param p Output array to store the counts of each bit pair.
 */
template <typename T>
void bit_pair_count_bits(bits<T> a, bits<T> b, int j, int *p) {
	int bit_a, bit_b;

	if (std::is_same_v<T, float>) {
		uint32_t mask = 0x00000001;
		mask <<= j;
		bit_a = (a.bit_set_32 & mask) ? 1 : 0;
		bit_b = (b.bit_set_32 & mask) ? 1 : 0;
	} else if (std::is_same_v<T, double>) {
		uint64_t mask = 0x0000000000000001;
		mask <<= j;
		bit_a = (a.bit_set_64 & mask) ? 1 : 0;
		bit_b = (b.bit_set_64 & mask) ? 1 : 0;
	}

	for (int i = 0; i < 4; i++) p[i] = 0;
	if (bit_a == 0 && bit_b == 0) {
		p[0] = 1;
	} else if (bit_a == 0 && bit_b == 1) {
		p[1] = 1;
	} else if (bit_a == 1 && bit_b == 0) {
		p[2] = 1;
	} else if (bit_a == 1 && bit_b == 1) {
		p[3] = 1;
	}
}

/**
 * @brief Counts the occurrences of bit pairs (00, 01, 10, 11) between two arrays.
 * @param A Pointer to the first array.	
 * @param B Pointer to the second array.
 * @param n_elem The number of elements in the arrays A and B.
 * @param pair_counts Pointer to an array to store the counts of each bit pair.
 */
template <typename T>
void bit_pair_count_bits(T *A, T *B, size_t n_elem, size_t *pair_counts) {
	bit_pair_count_calls += 1;
	const int n_bits = sizeof(T) * CHAR_BIT;
	int start_bit, end_bit;
	get_start_end_bits(n_bits, &start_bit, &end_bit);
	int total_bits = end_bit - start_bit;
	#pragma omp parallel for reduction(+:pair_counts[:total_bits * n_joint_counts])
	for (size_t i = 0; i < n_elem; i++) {
		bits<T> a(A[i]);
		bits<T> b(B[i]);
		int pair_count[4];
		for (int j = start_bit; j < end_bit; i++) {
			bit_pair_count_bits(a, b, j, pair_count);
			pair_counts[(j - start_bit) * n_joint_counts] += pair_count[0];
			pair_counts[(j - start_bit) * n_joint_counts + 1] += pair_count[1];
			pair_counts[(j - start_bit) * n_joint_counts + 2] += pair_count[2];
			pair_counts[(j - start_bit) * n_joint_counts + 3] += pair_count[3];
		}
	}
}

/**
 * @brief Helper function to compute pair counts at a particular position in two numbers represented as bitsets.
 * @param a First number represented as a bitset.
 * @param b Second number represented as a bitset.
 * @param j Bit position to compare.
 * @param p Output array to store the counts of each bit pair.
 */
template <int n_bits>
void bit_pair_count(std::bitset<n_bits> a, std::bitset<n_bits> b, int j, int *p) {
	for (int i = 0; i < 4; i++) p[i] = 0;
	if (a[j] == 0 && b[j] == 0) {
		p[0] = 1;
	} else if (a[j] == 0 && b[j] == 1) {
		p[1] = 1;
	} else if (a[j] == 1 && b[j] == 0) {
		p[2] = 1;
	} else if (a[j] == 1 && b[j] == 1) {
		p[3] = 1;
	}
}

/**
 * @brief Counts the occurrences of bit pairs (00, 01, 10, 11) between two arrays.
 *
 * @tparam T The type of the elements in the input arrays.
 * @param A Pointer to the first array.
 * @param B Pointer to the second array.
 * @param n_elem The number of elements in the arrays A and B.
 * @param pair_counts Pointer to an array to store the counts of each bit pair.
 */
template <typename T>
void bit_pair_count(T *A, T *B, size_t n_elem, size_t *pair_counts) {
	bit_pair_count_calls += 1;
	const int n_bits = sizeof(T) * CHAR_BIT;
	int start_bit, end_bit;
	get_start_end_bits(n_bits, &start_bit, &end_bit);
	int total_bits = end_bit - start_bit;
	#pragma omp parallel for reduction(+:pair_counts[:total_bits * n_joint_counts])
	for (size_t i = 0; i < n_elem; i++) {
		auto a = reinterpret_cast<typename HELP<T>::type *>(&A[i]);
		auto b = reinterpret_cast<typename HELP<T>::type *>(&B[i]);
		int pair_count[4];
		for (int j = start_bit; j < end_bit; j++) {
			bit_pair_count<sizeof(T) * CHAR_BIT>(*a, *b, j, pair_count);
			pair_counts[(j - start_bit) * n_joint_counts] += pair_count[0];
			pair_counts[(j - start_bit) * n_joint_counts + 1] += pair_count[1];
			pair_counts[(j - start_bit) * n_joint_counts + 2] += pair_count[2];
			pair_counts[(j - start_bit) * n_joint_counts + 3] += pair_count[3];
		}
	}
}
// End helper functions

/**
 * @brief Computes the entropy of each bit position in an array.
 *
 * @tparam T The type of the elements in the input array.
 * @param A Pointer to the array of numeric values.
 * @param n_elem The number of elements in the array A.
 * @param H Pointer to an array to store the entropy of each bit position.
 */
template <typename T>
void bit_count_entropy(T *A, size_t n_elem, double *H) {
	const int n_bits = sizeof(A[0]) * CHAR_BIT;
	int start_bit, end_bit;
	get_start_end_bits(n_bits, &start_bit, &end_bit);
	int total_bits = end_bit - start_bit;
	size_t c[total_bits];
	bit_count(A, n_elem, c);
	//bit_count_bits(A, n_elem, c); // bits class implementation

	for (int i = start_bit; i < end_bit; i++) {
		double p = static_cast<double>(c[i - start_bit])/static_cast<double>(n_elem);
		double p_array[2] = {p, 1 - p};
		H[i - start_bit] = entropy(p_array, 2);
	}
}

/** @brief shaves specified number of least significant bits from given bitset
 *
 *  @tparam n size of bitset
 *  @param b pointer to bitset to shave
 *  @param n number of bits to shave
 */
template <size_t nbits>
void shave_bitset(std::bitset<nbits> *b, int n) {
	std::bitset<nbits> mask = -1;	// all 1's
	mask <<= n;
	*b &= mask;
}

/**
 * @brief Clears the least significant bits of each element in an array.
 *
 * @tparam T The type of the elements in the input array.
 * @param a Pointer to the input array.
 * @param n_elem The number of elements in the input array.
 * @param n The number of least significant bits to clear.
 * @param shaved Pointer to the output array to store the modified elements.
 */
template <typename T>
void shave_template(T *a, size_t n_elem, int n, T *shaved) {
	for (size_t i = 0; i < n_elem; i++) {
		typename HELP<T>::type b;
		std::memcpy(&b, &a[i], sizeof(T));
		shave_bitset(&b, n);
		std::memcpy(&shaved[i], &b, sizeof(T));
	}
}

/** @brief replaces specified number of least significant bits from given bitset with 10000...
 *
 *  @tparam n size of bitset
 *  @param b pointer to bitset to shave
 *  @param n number of bits to shave
 */
template <size_t nbits>
void halfshave_bitset(std::bitset<nbits> *b, int n) {
	std::bitset<nbits> keep_mask = -1;	// all 1's
	keep_mask <<= n;
	*b &= keep_mask;
	std::bitset<nbits> round_up_mask = 0x80000000;
	round_up_mask >>= (nbits - n - 1);
	*b |= round_up_mask;
}

/**
 * @brief half shaves (i.e. replaces with 100...) the least significant bits of each element in an array.
 *
 * @tparam T The type of the elements in the input array.
 * @param a Pointer to the input array.
 * @param n_elem The number of elements in the input array.
 * @param n The number of least significant bits to clear.
 * @param halfshaved Pointer to the output array to store the modified elements.
 */
template <typename T>
void halfshave_template(T *a, size_t n_elem, int n, T *halfshaved) {
	for (size_t i = 0; i < n_elem; i++) {
		typename HELP<T>::type b;
		std::memcpy(&b, &a[i], sizeof(T));
		halfshave_bitset(&b, n);
		std::memcpy(&halfshaved[i], &b, sizeof(T));
	}
}

/**
 *  @brief Rounds the least significant bits of a bitset to the specified number of bits.
 *         Does this by creating an upper and lower mask to isolate the nth bit as a bitset
 *		   then adds it to the original bitset that had n-1 least significant bits cleared. This is a hacky way to do rounding since we can't do math operations with bitsets, but it works.
 *
 *  @tparam n size of bitset
 *  @param b pointer to bitset to round
 *  @param n number of bits to round to
 */
template <size_t nbits>
void round_bitset(std::bitset<nbits> *b, int n){
	int mantissa_bits = nbits - get_n_exponent_bits(nbits) - 1;
	std::bitset<nbits> mask = -1;	// all 1's
	std::bitset<nbits> lower_mask = mask << n;
	std::bitset<nbits> upper_mask = mask >> (nbits - n + 1);
	int keepbits = mantissa_bits - n;
	int shift = mantissa_bits - keepbits - 1;
	shift -= (keepbits == mantissa_bits) ? 1 : 0;

	if (nbits == 32) {
		std::uint32_t temp_b, upper_mask_uint;
		std::memcpy(&temp_b, b, sizeof(std::uint32_t));
		std::memcpy(&upper_mask_uint, &upper_mask, sizeof(std::uint32_t));
		temp_b += upper_mask_uint + ((temp_b >> shift) & 1);
		std::memcpy(b, &temp_b, sizeof(std::uint32_t));
		*b &= lower_mask;
	} else if (nbits == 64) {
		std::uint64_t temp_b, upper_mask_uint;
		std::memcpy(&temp_b, b, sizeof(std::uint64_t));
		std::memcpy(&upper_mask_uint, &upper_mask, sizeof(std::uint64_t));
		temp_b += upper_mask_uint + ((temp_b >> shift) & 1);
		std::memcpy(b, &temp_b, sizeof(std::uint64_t));
		*b &= lower_mask;
	} else {
		std::cerr << "Unsupported bit size: " << nbits << "\n";
	}
}

/**
 * @brief Performs bit rounding to the specified number of bits.
 *
 * @tparam T The type of the elements in the input array.
 * @param a Pointer to the input array.
 * @param n_elem The number of elements in the input array.
 * @param n The number of least significant bits to clear.
 * @param bit_rounding_result Pointer to the output array to store the modified elements.
 */
template <typename T>
void bit_rounding_template(T *a, size_t n_elem, int n, T *bit_rounding_result) {
	for (size_t i = 0; i < n_elem; i++) {
		typename HELP<T>::type b;
		std::memcpy(&b, &a[i], sizeof(T));
		round_bitset(&b, n);
		std::memcpy(&bit_rounding_result[i], &b, sizeof(T));
	}
}

/** @brief sets specified number of least significant bits from given bitset to 1
 *
 *  @tparam n size of bitset
 *  @param b pointer to bitset to set
 *  @param n number of bits to set
 */
template <size_t nbits>
void set_bitset(std::bitset<nbits> *b, int n) {
	std::bitset<nbits> mask = -1;	// all 1's
	mask <<= n;
	*b |= ~mask;
}

/**
 * @brief Clears the least significant bits of each element in an array.
 *
 * @tparam T The type of the elements in the input array.
 * @param a Pointer to the input array.
 * @param n_elem The number of elements in the input array.
 * @param n The number of least significant bits to clear.
 * @param set Pointer to the output array to store the modified elements.
 */
template <typename T>
void set_template(T *a, size_t n_elem, int n, T *set) {
	for (size_t i = 0; i < n_elem; i++) {
		typename HELP<T>::type b;
		std::memcpy(&b, &a[i], sizeof(T));
		set_bitset(&b, n);
		std::memcpy(&set[i], &b, sizeof(T));
	}
}

/**
 * @brief Clears the least significant bits of each element in an array.
 *
 * @tparam T The type of the elements in the input array.
 * @param a Pointer to the input array.
 * @param n_elem The number of elements in the input array.
 * @param n The number of least significant bits to clear.
 * @param groomed Pointer to the output array to store the modified elements.
 */
template <typename T>
void groom_template(T *a, size_t n_elem, int n, T *groomed) {
	for (size_t i = 0; i < n_elem; i++) {
		typename HELP<T>::type b;
		std::memcpy(&b, &a[i], sizeof(T));
		if (i % 2 == 0)
		  set_bitset(&b, n);
		else
		  shave_bitset(&b, n);
		std::memcpy(&groomed[i], &b, sizeof(T));
	}
}

/**
 * @brief Computes the mutual information for each bit position between two arrays.
 *
 * @tparam T The type of the elements in the input arrays.
 * @param A Pointer to the first input array.
 * @param B Pointer to the second input array.
 * @param n_elem The number of elements in the input arrays.
 * @param info Pointer to an array to store the mutual information for each bit position.
 */
template <typename T>
void mutual_information(T *A, T *B, size_t n_elem, double *info) {
	const int n_bits = sizeof(A[0]) * CHAR_BIT;
	int start_bit, end_bit;
	get_start_end_bits(n_bits, &start_bit, &end_bit);
	int total_bits = end_bit - start_bit;

	size_t pair_counts[n_joint_counts * total_bits];
	for (size_t i = 0; i < n_joint_counts * total_bits; i++) pair_counts[i] = 0;
	//bit_pair_count_bits(A, B, n_elem, pair_counts); // bits class implementation
	bit_pair_count(A, B, n_elem, pair_counts);

	for (int i = 0; i < total_bits; i++) {
		info[i] = 0;
		for (int r = 0; r < 2; r++) {
			double p_r = static_cast<double>(pair_counts[i * n_joint_counts + r * n_counts] + pair_counts[i * n_joint_counts + r * n_counts + 1]) / static_cast<double>(n_elem);
			for (int s = 0; s < 2; s++) {
				double p_s = static_cast<double>(pair_counts[i * n_joint_counts + s] + pair_counts[i * n_joint_counts + n_counts + s]) / static_cast<double>(n_elem);
				if (pair_counts[i * n_joint_counts + r * n_counts + s] > 0) {
					double p_rs = static_cast<double>(pair_counts[i * n_joint_counts + r * n_counts + s]) / static_cast<double>(n_elem);
					info[i] += p_rs * std::log2(p_rs/(p_r * p_s));
				}
			}
		}
	}
}

// XXX: doing 1D for now, depends how multidimensional arrays are stored. 
// If using mdspan things would be straightforward, however there is poor support (for now).
// If using pointer, need to do some fancy(ish) indexing.
/** 
 * @brief Compute real information stored in given array by evaluating mutual information of array with a shifted version of itself
 * 
 * @tparam T Type of elements in array
 * @param A Pointer to array
 * @param n_elem number of elements in A
 * @param info Pointer to array storing real information
 */
template <typename T>
void bitwise_real_information(T *A, size_t n_elem, double *info) {
	T *B = A + 1;
	mutual_information(A, B, n_elem - 1, info);
}

/**
 * @brief Computes amount of redundant information between two arrays
 *
 * @tparam The types of arrays A and B
 * @param A Pointer to first array 
 * @param B Pointer to first array 
 * @param n_elem Number of elements in each array 
 * @param R Pointer to array to store redundancy for each bit position 
 */
template <typename T>
void redundancy(T *A, T *B, size_t n_elem, double *R) {
	const int n_bits = sizeof(A[0]) * CHAR_BIT;
	int start_bit, end_bit;
	get_start_end_bits(n_bits, &start_bit, &end_bit);
	int total_bits = end_bit - start_bit;
	double *M = new double[total_bits];
	double *HA = new double[total_bits];
	double *HB = new double[total_bits];

	mutual_information(A, B, n_elem, M);
	bit_count_entropy(A, n_elem, HA);
	bit_count_entropy(B, n_elem, HB);

	for (int i = 0; i < total_bits; i++) {
		if (HA[i] + HB[i] > 0) {
			R[i] = 2*M[i]/(HA[i] + HB[i]);
		} else if (M[i] == HA[i] == HB[i] == 0) {
			R[i] = 1;  // The exponent bits are often identical between all elements, which mathematically works out to a redundancy of 1
		}
	}

}

/**
 *  @brief Computes real information preserved in array B with respect to A.
 *
 *  @tparam T Type of elements in arrays A and B
 *  @param A Pointer to original array
 *  @param B Pointer to shaved array
 *  @param n_elem Number of elements in each array
 *  @return Fraction of real information preserved in B with respect to A
 */
template <typename T>
double preserved_information_template(T *A, T *B, size_t n_elem) {
	const int n_bits = sizeof(A[0]) * CHAR_BIT;
	int start_bit, end_bit;
	get_start_end_bits(n_bits, &start_bit, &end_bit);
	int total_bits = end_bit - start_bit;
	double *R = new double[total_bits];
	double *I = new double[total_bits];

	redundancy(A, B, n_elem, R);
	bitwise_real_information(A, n_elem, I);

	double P = 0;
	double total_information = 0;
	for (int i = 0; i < total_bits; i++) {
		P += R[i] * I[i] ;
		total_information += I[i];
	}

    if (total_information == 0) {
		return 0.0; 
	}
	P /= total_information;
	return P;
}

/**
 * @brief Picks the number of bits to shave from elements in array A such that the preserved information is above a given tolerance.
 *
 * @tparam T The type of the elements in the input array.
 * @param A Pointer to the input array.
 * @param n_elem The number of elements in the input array.
 * @param tolerance The minimum fraction of preserved information required.
 * @param nbits_old The previous number of bits shaved (0 if starting fresh).
 * @return The number of bits to shave from each element in the array.
 */
template <typename T>
int pick_bits_to_shave_template(T *A, size_t n_elem, double tolerance, int nbits_old) {
  int max_bits = sizeof(T) * CHAR_BIT; 
  T s[n_elem];
  int start_bit;
  if constexpr (std::is_same_v<T, float>) { 
    start_bit = max_bits - n_exponent_bits_float - 1;
  } else if constexpr (std::is_same_v<T, double>) {
    start_bit = max_bits - n_exponent_bits_double - 1;
  }
  if (nbits_old == 0) {
    /* Start with shaving more bits rather than less */
    for (int i = start_bit; i >= 0; i--) {
    	shave_template<T>(A, n_elem, i, s);
    	//shave_bits<T>(A, n_elem, i, s); // bits class implementation
    	auto info = preserved_information_template<T>(A, s, n_elem);
    	if (info > tolerance) {
    		return i;
    	}
    }
  } else {
    int nbits = nbits_old;
    	
    /* Check if we should retain more bits first */
    shave_template<T>(A, n_elem, nbits, s);
    //shave_bits<T>(A, n_elem, nbits, s); // bits class implementation
    auto info = preserved_information_template<T>(A, s, n_elem);
    if (info <= tolerance) {
      for (int i = nbits - 1; i >= 0; i--) {
      	shave_template<T>(A, n_elem, i, s);
      	//shave_bits<T>(A, n_elem, i, s); // bits class implementation
      	auto info = preserved_information_template<T>(A, s, n_elem);
      	if (info > tolerance) {
      		return i;
      	}
      }
    } else {
  	/* Maybe we aren't throwing out enough. Increase number of bits thrown out */
  	for (int i = nbits + 1; i < max_bits; i++) {
      	shave_template<T>(A, n_elem, i, s);
      	//shave_bits<T>(A, n_elem, i, s);
      	auto info = preserved_information_template<T>(A, s, n_elem);
      	if (info <= tolerance) {
      		return i-1; // XXX: Not sure if this should be i-1 or just i...
      	}
  	}
    }

  }

  return 0;
}

/**
 * @brief Gets the maximum number of bits that can be shaved for a floating point type.
 *
 * @tparam T The floating point type (float or double).
 * @return The maximum number of least significant bits that can be cleared while
 *         preserving the exponent bits and sign bit.
 */
template <typename T>
int get_max_shave_bits() {
  int max_bits = sizeof(T) * CHAR_BIT; 
  if constexpr (std::is_same_v<T, float>) { 
    return max_bits - n_exponent_bits_float - 1;
  } else if constexpr (std::is_same_v<T, double>) {
    return max_bits - n_exponent_bits_double - 1;
  }
  return 0;
}

/**
 * @brief Binary search to find optimal number of bits to shave while preserving information.
 *
 * @tparam T The floating point type (float or double).
 * @param A Pointer to input array.
 * @param n_elem Number of elements in array.
 * @param tolerance Minimum information preservation threshold.
 * @param start_bit Initial number of bits to try shaving.
 * @param high_bit Upper bound for binary search.
 * @param low_bit Lower bound for binary search.
 * @return Number of bits that can be shaved while maintaining information above tolerance.
 */
template <typename T>
int binary_search(T *A, size_t n_elem, double tolerance, int start_bit, int high_bit, int low_bit) {
	T s_high[n_elem];
	T s_low[n_elem];
    shave_template<T>(A, n_elem, start_bit, s_high);
    shave_template<T>(A, n_elem, start_bit-1, s_low);
    //shave_bits<T>(A, n_elem, start_bit, s_high);  // bits class implementation
    //shave_bits<T>(A, n_elem, start_bit-1, s_low); // bits class implementation
    //bit_rounding_template<T>(A, n_elem, start_bit, s_high);
    //bit_rounding_template<T>(A, n_elem, start_bit-1, s_low);
   	auto info_high = preserved_information_template<T>(A, s_high, n_elem);
   	auto info_low = preserved_information_template<T>(A, s_low, n_elem);
	// Note: info_low is always >= info_high since we are shaving fewer bits
   	if (info_high < tolerance && info_low >= tolerance) {
   		return start_bit-1;
   	} else if (info_low < tolerance) {
		int new_start_bit = (start_bit - low_bit)/2 + low_bit;
		if (new_start_bit == 0) return new_start_bit;
		return binary_search<T>(A, n_elem, tolerance, new_start_bit, start_bit, low_bit);
	} else if (info_high >= tolerance) {
		int new_start_bit = (high_bit - start_bit)/2 + start_bit;
		if (new_start_bit == start_bit) return new_start_bit;
		return binary_search<T>(A, n_elem, tolerance, new_start_bit, high_bit, start_bit);
	} else {
		std::cout << "binary search error start bit " << start_bit << " old bits shaved " << start_bit << " info high " << info_high << " info low " << info_low << " tolerance " << tolerance << "\n";
		return 0;
	}
}

/**
 * @brief Picks the number of bits to shave from elements in array A such that the preserved information is above a given tolerance.
 *
 * @tparam T The type of the elements in the input array.
 * @param A Pointer to the input array.
 * @param n_elem The number of elements in the input array.
 * @param tolerance The minimum fraction of preserved information required.
 * @param nbits_old The previous number of bits shaved (0 if starting fresh).
 * @return The number of bits to shave from each element in the array.
 */
template <typename T>
int pick_bits_to_shave_binary_search_template(T *A, size_t n_elem, double tolerance, int nbits_old) {
  int max_shave_bits = get_max_shave_bits<T>();

  bit_pair_count_calls = 0;
  if (nbits_old == 0) {
	return binary_search<T>(A, n_elem, tolerance, max_shave_bits/2, max_shave_bits, 0);
  } else {
	return binary_search<T>(A, n_elem, tolerance, nbits_old, max_shave_bits, 0);
  }

  return 0;
}

#endif