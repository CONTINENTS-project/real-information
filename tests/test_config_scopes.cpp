#include "tests.h"
#include "bit_info.h"

#include <sstream>

yaml_config configuration(const std::string &text) {
    std::istringstream input(text);
    return yaml_config(input, "scopes.yaml");
}

template <typename Function>
void rejects(Function function, const std::string &detail) {
    bool rejected = false;
    try {
        function();
    } catch (const configuration_error &error) {
        rejected = true;
        const std::string message(error.what());
        test_assert(message.find("scopes.yaml") != std::string::npos &&
                    message.find(detail) != std::string::npos,
                    __FILE__, __LINE__, "Configuration error lacks context: " + message);
    }
    test_assert(rejected, __FILE__, __LINE__, "Expected configuration error: " + detail);
}

int main() {
    const auto config = configuration(
        "mantissa_only: false\n"
        "mask_zero: true\n"
        "use_real_info: true\n"
        "default_real_info: 0.99\n"
        "variables: # variable overrides\n"
        "  - var_name: PRECC\n"
        "    mantissa_only: true\n"
        "    mask_zero: false\n"
        "    real_info: 0.999\n"
        "\n"
        "  # Another variable; var_name need not be its first field.\n"
        "  - fixed_point: true\n"
        "    var_name: PM25_SRF\n"
        "fixed_point: false # global setting after the sequence\n");
    test_assert(!config.read<bool>("mantissa_only"), __FILE__, __LINE__, "Nested override leaked into globals");
    test_assert(!config.read<bool>("fixed_point"), __FILE__, __LINE__, "Global lookup did not resume after variables");
    test_assert(config.read<bool>("mantissa_only", "PRECC"), __FILE__, __LINE__, "Matching override was ignored");
    test_assert(!config.read<bool>("mask_zero", "PRECC"), __FILE__, __LINE__, "Explicit false override was ignored");
    test_assert(!config.read<bool>("fixed_point", "PRECC"), __FILE__, __LINE__, "Missing override did not inherit global");
    test_assert(config.read<bool>("fixed_point", "PM25_SRF"), __FILE__, __LINE__, "Later variable was not found");
    test_assert(config.read<bool>("mask_zero", "PM25_SRF"), __FILE__, __LINE__, "Inherited true was lost");
    test_assert(!config.read<bool>("mantissa_only", "UNLISTED"), __FILE__, __LINE__, "Unlisted variable did not inherit globals");
    test_assert(!config.read<bool>("mantissa_only", ""), __FILE__, __LINE__, "Empty variable name did not select globals");

    const auto missing = configuration(
        "variables:\n"
        "  - var_name: FIRST\n"
        "    mask_zero: true\n"
        "  - var_name: SECOND\n"
        "    fixed_point: true\n");
    rejects([&] { missing.read<bool>("mask_zero"); }, "global key 'mask_zero'");
    rejects([&] { missing.read<bool>("mask_zero", "SECOND"); }, "global key 'mask_zero'");
    rejects([&] { missing.read<bool>("mask_zero", "UNLISTED"); }, "global key 'mask_zero'");

    const auto quoted = configuration(
        "name: \"a#b:c\" # comment\r\n"
        "plain: a#b:c # comment\r\n"
        "apostrophe: O'Brien\r\n"
        "escaped: \"a\\\"b\\\\c\\t\"\r\n"
        "variables:\r\n"
        "  - var_name: 'O''Brien#field' # comment\r\n"
        "    mask_zero: false\r\n"
        "mask_zero: true\r\n");
    test_assert(quoted.read<std::string>("name") == "a#b:c", __FILE__, __LINE__, "Quoted hash/colon changed");
    test_assert(quoted.read<std::string>("plain") == "a#b:c", __FILE__, __LINE__, "Plain hash/colon changed");
    test_assert(quoted.read<std::string>("apostrophe") == "O'Brien", __FILE__, __LINE__, "Plain apostrophe changed");
    test_assert(quoted.read<std::string>("escaped") == "a\"b\\c\t", __FILE__, __LINE__, "Quoted escapes changed");
    test_assert(!quoted.read<bool>("mask_zero", "O'Brien#field"), __FILE__, __LINE__, "Quoted variable name did not match");
    test_assert(configuration("name: ''\n").read<std::string>("name").empty(),
                __FILE__, __LINE__, "Quoted empty string was rejected");

    const auto empty = configuration("variables: []\nmask_zero: true\n");
    test_assert(empty.read<bool>("mask_zero", "ANY"), __FILE__, __LINE__, "Empty sequence blocked globals");
    const auto next_line = configuration("variables:\n  -\n    var_name: A\n    mask_zero: true\n");
    test_assert(next_line.read<bool>("mask_zero", "A"), __FILE__, __LINE__, "Mapping after bare dash rejected");
    const auto indentless = configuration("variables:\n- var_name: A\n  mask_zero: true\n");
    test_assert(indentless.read<bool>("mask_zero", "A"), __FILE__, __LINE__, "Indentless sequence rejected");

    const std::vector<std::pair<std::string, std::string>> malformed = {
        {"mask_zero: false\nmask_zero: true\n", "duplicate key 'mask_zero'"},
        {"variables: []\nvariables: []\n", "duplicate key 'variables'"},
        {"variables:\n  - var_name: A\n  - var_name: A\n", "duplicate variable 'A'"},
        {"variables:\n  - var_name: A\n    mask_zero: true\n    mask_zero: false\n", "duplicate key 'mask_zero'"},
        {"variables:\n  - mask_zero: true\n", "var_name"},
        {"variables:\n  - var_name: ''\n", "nonempty var_name"},
        {"variables: {}\n", "variables must be a sequence"},
        {"variables: null\n", "variables must be a sequence"},
        {"variables:\n", "use []"},
        {"variables:\nmask_zero: true\n", "use []"},
        {"variables:\n  var_name: A\n", "unexpected indentation"},
        {"variables:\n  - A\n", "key: value"},
        {"- var_name: A\n", "outside variables"},
        {"variables:\n  - var_name: A\n   - var_name: B\n", "sequence indentation"},
        {"variables:\n  - var_name: A\n      mask_zero: true\n", "field indentation"},
        {"variables:\n  - var_name: A\n    nested:\n      mask_zero: true\n", "expected a scalar"},
        {"nested:\n  mask_zero: true\n", "expected a scalar"},
        {"  mask_zero: true\n", "unexpected indentation"},
        {"\tmask_zero: true\n", "tabs"},
        {"variables:\n  -\tvar_name: A\n", "tabs"},
        {"mask_zero:true\n", "key: value"},
        {"mask_zero\n", "key: value"},
        {"name: 'unfinished\n", "unterminated"},
        {"name: \"unfinished\\\n", "unfinished quoted escape"},
        {"name: \"unsupported\\q\"\n", "unsupported quoted escape"},
        {"name: 'closed' trailing\n", "unexpected text"},
        {"name: 'closed'#comment\n", "expected whitespace"},
        {"name: a: b\n", "nested mappings"},
        {"name: [a, b]\n", "unsupported YAML"},
        {"name: ]\n", "unsupported YAML"},
        {"name: &anchor true\n", "unsupported YAML"},
        {"mask_zero: 'true'\n", "unquoted scalar"},
        {"mask_zero: ture\n", "mask_zero"},
        {"variables:\n  - var_name: UNUSED\n    fixed_point: ture\n", "variables[UNUSED].fixed_point"},
    };
    for (const auto &[text, detail] : malformed) {
        rejects([&] { configuration(text); }, detail);
    }
    return 0;
}
