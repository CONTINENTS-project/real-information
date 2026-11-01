#include <bit_info.h>

double entropy_from_second_file() {
    double probabilities[] = {0.5, 0.5};
    return entropy(probabilities, 2);
}
