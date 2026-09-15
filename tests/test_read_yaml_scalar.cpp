#include "tests.h"
#include "bit_info.h"

int main() {
    // write path to test yaml
    const char *cfg = "/work/d446/d446/ab_continents/cesm_workflow/real-information/tests/test_yaml_override.yaml";

    // Global values
    bool mantissa_global = false;
    bool mask_zero_global = true;
    bool fixed_point_global = true;
    int fixed_point_scale_bits_global = 10;

    read_yaml_scalar(cfg, "mantissa_only", &mantissa_global);
    read_yaml_scalar(cfg, "mask_zero", &mask_zero_global);
    read_yaml_scalar(cfg, "fixed_point", &fixed_point_global);
    read_yaml_scalar(cfg, "fixed_point_scale_bits", &fixed_point_scale_bits_global);

    test_assert(mantissa_global == false, __FILE__, __LINE__, "mantissa_only global mismatch");
    test_assert(mask_zero_global == true, __FILE__, __LINE__, "mask_zero global mismatch");
    test_assert(fixed_point_global == true, __FILE__, __LINE__, "fixed_point global mismatch");
    test_assert(fixed_point_scale_bits_global == 10, __FILE__, __LINE__, "fixed_point_scale_bits global mismatch");

    // Per-variable override: PRECC should override mantissa_only -> true and fixed_point_scale_bits -> 5
    bool mantissa_precc = false;
    int fixed_point_scale_precc = 10;
    read_yaml_scalar(cfg, "mantissa_only", &mantissa_precc, "PRECC");
    read_yaml_scalar(cfg, "fixed_point_scale_bits", &fixed_point_scale_precc, "PRECC");

    test_assert(mantissa_precc == true, __FILE__, __LINE__, "mantissa_only PRECC override failed");
    test_assert(fixed_point_scale_precc == 5, __FILE__, __LINE__, "fixed_point_scale_bits PRECC override failed");

    // Per-variable override: PM25_SRF should override mantissa_only -> true and fixed_point -> false
    bool mantissa_pm25 = false;
    bool fixed_point_pm25 = true;
    read_yaml_scalar(cfg, "mantissa_only", &mantissa_pm25, "PM25_SRF");
    read_yaml_scalar(cfg, "fixed_point", &fixed_point_pm25, "PM25_SRF");

    test_assert(mantissa_pm25 == true, __FILE__, __LINE__, "mantissa_only PM25_SRF override failed");
    test_assert(fixed_point_pm25 == false, __FILE__, __LINE__, "fixed_point PM25_SRF override failed");

    return 0;
}
