#include <bit_info.h>

#ifndef _OPENMP
#error "The installed target must propagate OpenMP compilation requirements"
#endif
static_assert(__cplusplus >= 202002L, "The installed target must require C++20");

double entropy_from_second_file();

int main(int argc, char **argv) {
    if (argc != 2) return 1;
    double values[] = {1.0, 1.5, 2.0, 2.5};
    real_info_params<double> params(argv[1], values, 4);
    std::array<double, 64> information{};
    mutual_information(params, values, values, 4, information.data());
    double total = 0;
    for (double value : information) total += value;
    if (!(total > 0) || entropy_from_second_file() != 1.0) return 1;
    // Also use an exported C symbol, ensuring the shared library is linked/loaded.
    return real_info_last_fixed_point_error() == nullptr ? 0 : 1;
}
