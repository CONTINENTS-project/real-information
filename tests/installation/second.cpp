#include <my_rounding.h>

double entropy_from_second_file() {
    double probabilities[] = {0.5, 0.5};
    return entropy(probabilities, 2);
}
