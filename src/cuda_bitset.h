#ifndef CUDA_BITSET_H
#define CUDA_BITSET_H

#include <stdio.h>
#include <cuda_runtime.h>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
  if (code != cudaSuccess) 
  {
    fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
    if (abort) exit(code);
  }
}

template <typename T>
class bits {
  public:
    uint bit_set;
    long long bit_set_ll;

	  __device__ bits() = default;

    __device__ bits(T x) {
      if constexpr (std::is_same_v<T, float>) {
        bit_set = __float_as_uint(x);
      } else if constexpr (std::is_same_v<T, double>) {
        bit_set_ll = __double_as_longlong(x);
      } else {
        // Abort
        static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>, "Unsupported bit size");
      }
    };

    __device__ T convert_back() {
      if constexpr (std::is_same_v<T, float>) {
        return __uint_as_float(bit_set);
      } else if constexpr (std::is_same_v<T, double>) {
        return __longlong_as_double(bit_set_ll);
      } else {
        // Abort
        static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>, "Unsupported bit size");
      }
    };

    __device__ int operator[](int n) {
      if constexpr (std::is_same_v<T, float>) {
        return (bit_set >> (sizeof(T) * CHAR_BIT - n - 1)) & 1; // alternative way to get the nth bit.
      } else if constexpr (std::is_same_v<T, double>) {
        return (bit_set_ll >> (sizeof(T) * CHAR_BIT - n - 1)) & 1; // alternative way to get the nth bit.
      } else {
        // Abort
        static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>, "Unsupported bit size");
      }
    } ;   

    __device__ void operator<<=(int n){
      if constexpr (std::is_same_v<T, float>) {
        bit_set <<= n;
      } else if constexpr (std::is_same_v<T, double>) {
        bit_set_ll <<= n;
      } else {
        // Abort
        static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>, "Unsupported bit size");
      }
    }
    
    __device__ void operator>>=(int n){
      if constexpr (std::is_same_v<T, float>) {
        bit_set >>= n;
      } else if constexpr (std::is_same_v<T, double>) {
        bit_set_ll >>= n;
      } else {
        // Abort
        static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>, "Unsupported bit size");
      }
    }
    
    __device__ void operator&=(bits<T> b){
      if constexpr (std::is_same_v<T, float>) {
        bit_set &= b.bit_set;
      } else if constexpr (std::is_same_v<T, double>) {
        bit_set_ll &= b.bit_set_ll;
      } else {
        // Abort
        static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>, "Unsupported bit size");
      }
    }
};

#endif