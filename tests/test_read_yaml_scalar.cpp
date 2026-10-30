#include "tests.h"
#include "bit_info.h"

int main() {
    const std::string cfg = test_data_path("test_yaml_override.yaml");

    // Global values
    bool mantissa_global = true;
    bool mask_zero_global = false;
    bool fixed_point_global = false;

    read_yaml_scalar(cfg.c_str(), "mantissa_only", &mantissa_global);
    read_yaml_scalar(cfg.c_str(), "mask_zero", &mask_zero_global);
    read_yaml_scalar(cfg.c_str(), "fixed_point", &fixed_point_global);

    test_assert(mantissa_global == false, __FILE__, __LINE__, "mantissa_only global mismatch");
    test_assert(mask_zero_global == true, __FILE__, __LINE__, "mask_zero global mismatch");
    test_assert(fixed_point_global == true, __FILE__, __LINE__, "fixed_point global mismatch");

    // PRECC overrides mantissa_only -> true and explicitly overrides mask_zero -> false.
    bool mantissa_precc = false;
    bool mask_zero_precc = true;
    read_yaml_scalar(cfg.c_str(), "mantissa_only", &mantissa_precc, "PRECC");
    read_yaml_scalar(cfg.c_str(), "mask_zero", &mask_zero_precc, "PRECC");

    test_assert(mantissa_precc == true, __FILE__, __LINE__, "mantissa_only PRECC override failed");
    test_assert(!mask_zero_precc, __FILE__, __LINE__, "mask_zero PRECC override failed");

    // Per-variable override: PM25_SRF should override mantissa_only -> true and fixed_point -> false
    bool mantissa_pm25 = false;
    bool fixed_point_pm25 = true;
    read_yaml_scalar(cfg.c_str(), "mantissa_only", &mantissa_pm25, "PM25_SRF");
    read_yaml_scalar(cfg.c_str(), "fixed_point", &fixed_point_pm25, "PM25_SRF");

    test_assert(mantissa_pm25 == true, __FILE__, __LINE__, "mantissa_only PM25_SRF override failed");
    test_assert(fixed_point_pm25 == false, __FILE__, __LINE__, "fixed_point PM25_SRF override failed");

    return 0;
}
