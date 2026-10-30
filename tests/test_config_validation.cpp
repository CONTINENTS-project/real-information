#include "tests.h"
#include "bit_info.h"

#include <sstream>

template <typename Function>
void expect_configuration_error(Function function, const std::string &detail) {
    bool rejected = false;
    try {
        function();
    } catch (const configuration_error &error) {
        rejected = true;
        test_assert(std::string(error.what()).find(detail) != std::string::npos,
                    __FILE__, __LINE__, "Configuration error lacks context: " + std::string(error.what()));
    }
    test_assert(rejected, __FILE__, __LINE__, "Expected configuration error: " + detail);
}

template <typename T>
T read_value(const std::string &text) {
    std::istringstream input("value: " + text + "\n");
    const yaml_config config(input, "conversion.yaml");
    T value{};
    read_yaml_scalar(config, "value", &value);
    return value;
}

template <typename T>
void reject_value(const std::string &text) {
    expect_configuration_error([&] { read_value<T>(text); }, "conversion.yaml:1: value");
}

int main() {
    for (const std::string value : {"true", "True", "TRUE"}) {
        test_assert(read_value<bool>(value), __FILE__, __LINE__, "Valid true rejected");
    }
    for (const std::string value : {"false", "False", "FALSE"}) {
        test_assert(!read_value<bool>(value), __FILE__, __LINE__, "Valid false rejected");
    }
    for (const std::string value : {"", "ture", "yes", "0", "false trailing"}) reject_value<bool>(value);

    test_assert(read_value<int>("+10") == 10, __FILE__, __LINE__, "Positive integer rejected");
    test_assert(read_value<int>("-10") == -10, __FILE__, __LINE__, "Negative integer rejected");
    test_assert(read_value<int8_t>("-128") == -128, __FILE__, __LINE__, "Lower integer endpoint rejected");
    test_assert(read_value<uint8_t>("255") == 255, __FILE__, __LINE__, "Upper integer endpoint rejected");
    for (const std::string value : {"", "10junk", "1.5", "+-1", "++1", "99999999999999999999999"}) {
        reject_value<int>(value);
    }
    reject_value<int8_t>("128");
    reject_value<int8_t>("-129");
    reject_value<uint8_t>("256");
    reject_value<unsigned>("-1");

    test_assert(read_value<float>("1.25") == 1.25f, __FILE__, __LINE__, "Finite float rejected");
    test_assert(read_value<double>("-2e3") == -2000.0, __FILE__, __LINE__, "Scientific notation rejected");
    test_assert(std::signbit(read_value<double>("-0.0")), __FILE__, __LINE__, "Signed zero lost");
    for (const std::string value : {"", "1.5junk", "nan", "inf", "-inf", "1e999", "1e-999"}) {
        reject_value<double>(value);
    }
    reject_value<float>("1e39");
    test_assert(read_value<std::string>("PRECC") == "PRECC", __FILE__, __LINE__, "String scalar rejected");

    expect_configuration_error([] { yaml_config config(nullptr); }, "path");
    expect_configuration_error([] { yaml_config config(""); }, "path");
    expect_configuration_error([] {
        yaml_config config(test_data_path("nonexistent_config_validation_fixture.yaml").c_str());
    }, "could not open");
    expect_configuration_error([] {
        std::istringstream input;
        input.setstate(std::ios::badbit);
        yaml_config config(input, "broken-stream.yaml");
    }, "broken-stream.yaml");

    std::istringstream input("value: invalid\n");
    const yaml_config config(input, "unchanged.yaml");
    int output = 123;
    expect_configuration_error([&] { read_yaml_scalar(config, "value", &output); }, "value");
    test_assert(output == 123, __FILE__, __LINE__, "Failed conversion changed its output");
    expect_configuration_error([&] { read_yaml_scalar(config, "missing", &output); }, "missing");
    expect_configuration_error([&] { read_yaml_scalar<int>(config, "value", nullptr); }, "null");

    std::string path = test_data_path("real_info.yaml");
    real_info_params<float> no_data(path.data());
    float data[] = {1.0f, 2.0f};
    real_info_params<float> with_data(path.data(), data, 2);
    test_assert(!no_data.mantissa_only && !no_data.mask_zero && !no_data.fixed_point,
                __FILE__, __LINE__, "Data-free constructor loaded wrong flags");
    test_assert(!with_data.mantissa_only && !with_data.mask_zero && !with_data.fixed_point,
                __FILE__, __LINE__, "Data constructor loaded wrong flags");
    test_assert(no_data.fixed_range.constant && no_data.fixed_range.width == 0,
                __FILE__, __LINE__, "Inactive fixed-point metadata is not initialized");
    path = test_data_path("missing_config_key.yaml");
    expect_configuration_error([&] { real_info_params<float> params(path.data()); }, "fixed_point");
    expect_configuration_error([&] { real_info_params<float> params(path.data(), data, 2); }, "fixed_point");
    return 0;
}
