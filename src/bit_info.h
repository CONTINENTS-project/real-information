#ifndef BIT_INFO_H
#define BIT_INFO_H

#include <algorithm>
#include <array>
#include <bit>
#include <charconv>
#include <concepts>
#include <climits>
#include <limits>
#include <iostream>
#include <fstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include <unordered_map>
#include <type_traits>
#include <cmath>
#include <cassert>
#include <cstring>
#include <cstdint>
#include <stdexcept>

// Set EPS
static const float EPS = 1e-16;

inline std::string yaml_trim(const std::string &s) {
	size_t start = s.find_first_not_of(" \t\r\n");
	if (start == std::string::npos) return "";
	size_t end = s.find_last_not_of(" \t\r\n");
	return s.substr(start, end - start + 1);
}

class configuration_error : public std::invalid_argument {
public:
	using std::invalid_argument::invalid_argument;
};

// Convert the entire scalar directly into its destination type, without narrowing.
template <typename T>
T parse_yaml_scalar(const std::string &value, const std::string &location) {
	if constexpr (std::is_same_v<T, std::string>) {
		return value;
	} else if constexpr (std::is_same_v<T, bool>) {
		if (value == "true" || value == "True" || value == "TRUE") return true;
		if (value == "false" || value == "False" || value == "FALSE") return false;
	} else if constexpr (std::is_integral_v<T> || std::is_floating_point_v<T>) {
		std::string_view text(value);
		if (!text.empty() && text.front() == '+') {
			text.remove_prefix(1);
			if (text.empty() || text.front() == '+' || text.front() == '-') {
				throw configuration_error(location + ": invalid numeric scalar '" + value + "'");
			}
		}
		if (!text.empty()) {
			T parsed{};
			const auto result = std::from_chars(text.data(), text.data() + text.size(), parsed);
			if (result.ec == std::errc{} && result.ptr == text.data() + text.size()) {
				if constexpr (std::is_integral_v<T>) return parsed;
				else if (std::isfinite(parsed)) return parsed;
			}
		}
	} else {
		static_assert(!sizeof(T), "Unsupported type for read_yaml_scalar");
	}
	throw configuration_error(location + ": invalid scalar '" + value + "'");
}

// Dependency-free subset: top-level scalars and a variables sequence of scalar maps.
// Unsupported YAML structures are rejected rather than searched as flat text.
class yaml_config {
    struct scalar {
        std::string value;
        size_t line;
        bool quoted = false;
    };
    using scalar_map = std::unordered_map<std::string, scalar>;
    scalar_map globals;
    std::unordered_map<std::string, scalar_map> variables;

    std::string location(size_t line) const {
        return filename + ":" + std::to_string(line);
    }

    [[noreturn]] void fail(size_t line, const std::string &message) const {
        throw configuration_error(location(line) + ": " + message);
    }

    static bool space(char character) {
        return character == ' ' || character == '\t';
    }

    scalar parse_scalar(const std::string &text, size_t line) const {
        scalar result{"", line};
        if (text.empty() || text.front() == '#') return result;
        if (text.front() == '\'' || text.front() == '"') {
            result.quoted = true;
            const char quote = text.front();
            for (size_t i = 1; i < text.size(); ++i) {
                char character = text[i];
                if (character == quote) {
                    if (quote == '\'' && i + 1 < text.size() && text[i + 1] == '\'') {
                        result.value += '\'';
                        ++i;
                        continue;
                    }
                    if (i + 1 < text.size() && !space(text[i + 1])) {
                        fail(line, "expected whitespace after quoted scalar");
                    }
                    const auto trailing = yaml_trim(text.substr(i + 1));
                    if (!trailing.empty() && trailing.front() != '#') {
                        fail(line, "unexpected text after quoted scalar");
                    }
                    return result;
                }
                if (quote == '"' && character == '\\') {
                    if (++i == text.size()) fail(line, "unfinished quoted escape");
                    switch (text[i]) {
                        case '"': character = '"'; break;
                        case '\\': character = '\\'; break;
                        case 'n': character = '\n'; break;
                        case 'r': character = '\r'; break;
                        case 't': character = '\t'; break;
                        default: fail(line, "unsupported quoted escape");
                    }
                }
                result.value += character;
            }
            fail(line, "unterminated quoted scalar");
        }
        size_t comment = std::string::npos;
        for (size_t i = 0; i < text.size(); ++i) {
            if (text[i] == '#' && (i == 0 || space(text[i - 1]))) {
                comment = i;
                break;
            }
        }
        result.value = yaml_trim(text.substr(0, comment));
        return result;
    }

    std::pair<std::string, scalar> parse_entry(const std::string &text, size_t line) const {
        const size_t colon = text.find(':');
        if (colon == std::string::npos || (colon + 1 < text.size() && !space(text[colon + 1]))) {
            fail(line, "expected 'key: value'");
        }
        const std::string key = yaml_trim(text.substr(0, colon));
        if (key.empty() || key.find_first_not_of(
                "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-") != std::string::npos) {
            fail(line, "expected an unquoted scalar key");
        }
        return {key, parse_scalar(yaml_trim(text.substr(colon + 1)), line)};
    }

    template <typename T>
    T convert(const scalar &value, const std::string &key, const std::string &scope) const {
        const std::string context = location(value.line) + ": " + scope + key;
        if constexpr (!std::is_same_v<T, std::string>) {
            if (value.quoted) throw configuration_error(context + ": expected an unquoted scalar");
        }
        return parse_yaml_scalar<T>(value.value, context);
    }

    void insert(scalar_map &fields, const std::string &key, scalar value) const {
        if (key == "fixed_point_scale_bits") {
            fail(value.line, "fixed_point_scale_bits has been removed; fixed-point uses a fixed 63-bit range mapping");
        }
        if (!value.quoted) {
            if (value.value.empty()) fail(value.line, key + ": expected a scalar value");
            if (value.value.find_first_of("[{}]|>&*!%@`") == 0 || value.value == "~" ||
                value.value == "null" || value.value == "Null" || value.value == "NULL" ||
                value.value == "-" || value.value.compare(0, 2, "- ") == 0) {
                fail(value.line, key + ": unsupported YAML scalar or structure");
            }
            for (size_t i = 0; i < value.value.size(); ++i) {
                if (value.value[i] == ':' && (i + 1 == value.value.size() || space(value.value[i + 1]))) {
                    fail(value.line, key + ": nested mappings require an unsupported YAML structure");
                }
            }
        }
        const auto inserted = fields.emplace(key, value);
        if (!inserted.second) fail(value.line, "duplicate key '" + key + "'");
    }

    void validate_flags(const scalar_map &fields, const std::string &scope) const {
        for (const char *key : {"mantissa_only", "mask_zero", "fixed_point"}) {
            const auto found = fields.find(key);
            if (found != fields.end()) convert<bool>(found->second, key, scope);
        }
    }

    void load(std::istream &input) {
        const size_t unset = std::string::npos;
        bool in_variables = false;
        bool seen_variables = false;
        bool have_item = false;
        size_t section_line = 0;
        size_t item_line = 0;
        size_t item_indent = unset;
        size_t field_indent = unset;
        scalar_map item;
        auto finish_item = [&] {
            if (!have_item) return;
            const auto name = item.find("var_name");
            if (name == item.end() || yaml_trim(name->second.value).empty()) {
                fail(item_line, "each variables entry needs a nonempty var_name");
            }
            const auto variable = name->second.value;
            validate_flags(item, "variables[" + variable + "].");
            if (!variables.emplace(variable, std::move(item)).second) {
                fail(item_line, "duplicate variable '" + variable + "'");
            }
            item.clear();
            have_item = false;
        };
        auto finish_section = [&] {
            if (in_variables && item_indent == unset) {
                fail(section_line, "variables must contain a sequence; use [] for an empty sequence");
            }
            finish_item();
            in_variables = false;
        };
        std::string line;
        size_t line_number = 0;
        while (std::getline(input, line)) {
            ++line_number;
            if (!line.empty() && line.back() == '\r') line.pop_back();
            const size_t indent = line.find_first_not_of(' ');
            if (indent == unset) continue;
            if (line[indent] == '\t') fail(line_number, "tabs are not allowed in indentation");
            const auto content = yaml_trim(line.substr(indent));
            if (content.empty() || content.front() == '#') continue;
            const bool sequence_item = content.front() == '-' &&
                (content.size() == 1 || space(content[1]));
            if (sequence_item) {
                if (!in_variables) fail(line_number, "sequence item outside variables");
                if (item_indent == unset) item_indent = indent;
                if (indent != item_indent) fail(line_number, "inconsistent variables sequence indentation");
                finish_item();
                have_item = true;
                item_line = line_number;
                field_indent = unset;
                const auto field_start = line.find_first_not_of(' ', indent + 1);
                if (field_start == unset || line[field_start] == '#') continue;
                if (line[field_start] == '\t') fail(line_number, "tabs are not allowed in indentation");
                field_indent = field_start;
                const auto entry = parse_entry(yaml_trim(line.substr(field_start)), line_number);
                insert(item, entry.first, entry.second);
            } else if (indent == 0) {
                finish_section();
                const auto entry = parse_entry(content, line_number);
                if (entry.first == "variables") {
                    if (seen_variables) fail(line_number, "duplicate key 'variables'");
                    seen_variables = true;
                    if (entry.second.quoted || (!entry.second.value.empty() && entry.second.value != "[]")) {
                        fail(line_number, "variables must be a sequence of scalar mappings");
                    }
                    in_variables = entry.second.value.empty();
                    section_line = line_number;
                    item_indent = unset;
                } else {
                    insert(globals, entry.first, entry.second);
                }
            } else {
                if (!in_variables || !have_item || indent <= item_indent) {
                    fail(line_number, "unexpected indentation outside a variables entry");
                }
                if (field_indent == unset) field_indent = indent;
                if (indent != field_indent) fail(line_number, "inconsistent variable field indentation");
                const auto entry = parse_entry(content, line_number);
                insert(item, entry.first, entry.second);
            }
        }
        if (input.bad() || (input.fail() && !input.eof())) {
            throw configuration_error(filename + ": failed to read configuration file");
        }
        finish_section();
        validate_flags(globals, "");
    }

public:
    std::string filename;

    explicit yaml_config(const char *path) {
        if (path == nullptr || *path == '\0') {
            throw configuration_error("Configuration path must be nonempty");
        }
        filename = path;
        std::ifstream file(path);
        if (!file.is_open()) throw configuration_error(filename + ": could not open configuration file");
        load(file);
    }

    yaml_config(std::istream &input, const std::string &source) : filename(source) {
        load(input);
    }

    template <typename T>
    T read(const std::string &key, const char *var_name = nullptr) const {
        if (var_name != nullptr && *var_name != '\0') {
            const auto variable = variables.find(var_name);
            if (variable != variables.end()) {
                const auto value = variable->second.find(key);
                if (value != variable->second.end()) {
                    return convert<T>(value->second, key, "variables[" + std::string(var_name) + "].");
                }
            }
        }
        const auto value = globals.find(key);
        if (value != globals.end()) return convert<T>(value->second, key, "");
        throw configuration_error(filename + ": required global key '" + key + "' not found");
    }
};

// Output remains unchanged on errors. Null/empty names select only global settings.
template <typename T>
void read_yaml_scalar(const yaml_config &config, const std::string &key, T *out, const char *var_name = nullptr) {
    if (out == nullptr) throw configuration_error(config.filename + ": null scalar output for " + key);
    *out = config.read<T>(key, var_name);
}

template <typename T>
void read_yaml_scalar(const char *config_file, const std::string &key, T *out, const char *var_name = nullptr) {
    const yaml_config config(config_file);
    read_yaml_scalar(config, key, out, var_name);
}

class fixed_point_error : public std::invalid_argument {
public:
	using std::invalid_argument::invalid_argument;
};

// Guarded C entry points print fixed-point validation errors and terminate
// the process. This accessor is retained for existing bindings and returns
// nullptr after successful calls.
extern "C" const char *real_info_last_fixed_point_error();

template <typename T>
struct fixed_point_range {
	T minimum{};
	T maximum{};
	double width = 0;
	bool constant = true;
};

template <typename T>
fixed_point_range<T> get_fixed_point_init_params(const T *x, size_t n_elem) {
	static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>);
	if (n_elem == 0 || x == nullptr) {
		throw fixed_point_error("Fixed-point initialization requires a nonempty array.");
	}
	for (size_t i = 0; i < n_elem; ++i) {
		if (!std::isfinite(x[i])) {
			throw fixed_point_error("Fixed-point input must contain only finite values.");
		}
	}
	const auto extrema = std::minmax_element(x, x + n_elem);
	fixed_point_range<T> range;
	range.minimum = *extrema.first;
	range.maximum = *extrema.second;
	range.constant = range.minimum == range.maximum;
	if (!range.constant) {
		// Subtract after promotion so a float range can exceed FLT_MAX safely.
		range.width = static_cast<double>(range.maximum) - static_cast<double>(range.minimum);
		if (!std::isfinite(range.width) || range.width <= 0) {
			throw fixed_point_error("Fixed-point range must have a finite positive width in double precision.");
		}
	}
	return range;
}

// Empty analysis inputs are invalid; singleton finite arrays are constant.
template <typename T>
bool is_constant_array(const T *x, size_t n_elem) {
	if (n_elem == 0 || x == nullptr) {
		throw fixed_point_error("Information analysis requires a nonempty array.");
	}
	return std::isfinite(x[0]) &&
		std::all_of(x + 1, x + n_elem, [first = x[0]](T value) { return value == first; });
}

template <typename T>
class real_info_params {
	private:
	    char *config_file = nullptr;

		void load_configuration() {
			const yaml_config config(config_file);
			read_yaml_scalar(config, "mantissa_only", &mantissa_only);
			read_yaml_scalar(config, "mask_zero", &mask_zero);
			read_yaml_scalar(config, "fixed_point", &fixed_point);
		}

	public:
	    bool mantissa_only = false;
		bool mask_zero = false;
		bool fixed_point = false;
		fixed_point_range<T> fixed_range{};
	
		real_info_params(char *config_file) {
			this->config_file = config_file;
			load_configuration();
			if (fixed_point) {
				throw fixed_point_error("Fixed-point initialization requires input data.");
			}
		}

		real_info_params(char *config_file, T *data, size_t n_elem) {
			this->config_file = config_file;
			load_configuration();
			if (fixed_point) {
				fixed_range = get_fixed_point_init_params(data, n_elem);
			}
		}
};

template <typename T>
bool copy_constant_fixed_point(const real_info_params<T> &params, const T *input,
	size_t n_elem, T *output) {
	if (!params.fixed_point || !params.fixed_range.constant) return false;
	if (n_elem == 0) return true;
	if (input == nullptr || output == nullptr) {
		throw fixed_point_error("Fixed-point conversion requires valid input and output arrays.");
	}
	for (size_t i = 0; i < n_elem; ++i) {
		if (!std::isfinite(input[i]) || input[i] != params.fixed_range.minimum) {
			throw fixed_point_error("Fixed-point value does not match the initialized constant range.");
		}
	}
	// Copy the original values to preserve the signs of zero as well.
	if (input != output) std::copy_n(input, n_elem, output);
	return true;
}

/*
  Fixed-point arithmetic class modelled on Alex Pinard's implementation.

  The endpoints map to 0 and max_code. max_code = 2^63 - 1024 is exactly
  representable as both double and int64_t, unlike double(INT64_MAX).
  Interior values are normalized in double precision and truncated.
  Constant ranges use code 0 and decode directly to their stored value.
*/
template <typename T>
class FP {
	static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>);
	static_assert(std::numeric_limits<double>::digits >= 53,
		"The fixed-point endpoint requires at least 53 bits of double precision.");
	fixed_point_range<T> range_;

	int64_t encode(double x) const {
		if (!std::isfinite(x) || x < static_cast<double>(range_.minimum) ||
			x > static_cast<double>(range_.maximum)) {
			throw fixed_point_error("Fixed-point value must be finite and inside the initialized range.");
		}
		if (range_.constant || x == static_cast<double>(range_.minimum)) return 0;
		if (x == static_cast<double>(range_.maximum)) return max_code;
		const double normalized = std::clamp(
			(x - static_cast<double>(range_.minimum)) / range_.width, 0.0, 1.0);
		return static_cast<int64_t>(normalized * static_cast<double>(max_code));
	}

	double normalized_value() const {
		if (value < 0 || value > max_code) {
			throw fixed_point_error("Fixed-point code is outside the initialized mapping.");
		}
		return static_cast<double>(value) / static_cast<double>(max_code);
	}

public:
	static constexpr int64_t max_code = std::numeric_limits<int64_t>::max() - 1023;
	int64_t value = 0;

	// Construct an encoded zero directly; do not try to encode floating-point 0.
	explicit FP(const fixed_point_range<T> &range) : range_(range) {}

	FP(const fixed_point_range<T> &range, T x) : range_(range) {
		value = encode(static_cast<double>(x));
	}

	T convert_back() const {
		if (range_.constant) return range_.minimum;
		const double normalized = normalized_value();
		if (value == 0) return range_.minimum;
		if (value == max_code) return range_.maximum;
		return static_cast<T>(std::lerp(static_cast<double>(range_.minimum),
			static_cast<double>(range_.maximum), normalized));
	}

	// Keep differences in double precision until the final output conversion.
	double convert_back_no_offset() const {
		if (range_.constant) return 0.0;
		return normalized_value() * range_.width;
	}

	FP operator+(const FP &other) const {
		FP result(range_);
		result.value = encode(static_cast<double>(convert_back()) + static_cast<double>(other.convert_back()));
		return result;
	}

	FP operator-(const FP &other) const {
		FP result(range_);
		result.value = encode(static_cast<double>(convert_back()) - static_cast<double>(other.convert_back()));
		return result;
	}
};

template <typename T>
void get_start_end_bits(const real_info_params<T> &params, int n_bits, int *start, int *end) {
	if (params.fixed_point) {
		*start = 0;
		*end = 63; // Assuming 64-bit fixed-point representation
	} else if (params.mantissa_only) {
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

namespace real_info_detail {

enum class bit_operation { shave, set, halfshave, round };

template <typename T>
using native_word = std::conditional_t<std::is_same_v<T, float>, uint32_t, uint64_t>;

template <typename T>
void validate_bit_count(int n, bool fixed_point = false) {
    static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>);
    const int maximum = fixed_point ? 63 : std::numeric_limits<T>::digits - 1;
    if (n < 0 || n > maximum) {
        throw fixed_point_error("Bit count must be between 0 and " + std::to_string(maximum));
    }
}

// Shared unsigned implementation; callers restrict native operations to the
// mantissa and fixed-point operations to the 63 nonnegative code bits.
template <std::unsigned_integral UInt>
UInt transform_word(UInt value, int n, bit_operation operation) {
    if (n < 0 || n >= std::numeric_limits<UInt>::digits) {
        throw fixed_point_error("Bit count is outside the unsigned word width");
    }
    if (n == 0) return value;
    const UInt step = UInt{1} << n;
    const UInt mask = step - UInt{1};
    const UInt cleared = value & ~mask;
    switch (operation) {
        case bit_operation::shave:
            return cleared;
        case bit_operation::set:
            return value | mask;
        case bit_operation::halfshave:
            return cleared | (UInt{1} << (n - 1));
        case bit_operation::round: {
            const UInt discarded = value & mask;
            const UInt halfway = step >> 1;
            // At a tie, the least-significant RETAINED bit decides parity.
            const bool increment = discarded > halfway ||
                (discarded == halfway && ((value >> n) & UInt{1}) != 0);
            if (!increment) return cleared;
            if (cleared > std::numeric_limits<UInt>::max() - step) {
                throw fixed_point_error("Rounded unsigned word would overflow");
            }
            return cleared + step;
        }
    }
    throw fixed_point_error("Unknown bit operation");
}

template <typename T>
native_word<T> transform_native(native_word<T> value, int n, bit_operation operation) {
    validate_bit_count<T>(n);
    using UInt = native_word<T>;
    constexpr int fraction_bits = std::numeric_limits<T>::digits - 1;
    constexpr UInt fraction_mask = (UInt{1} << fraction_bits) - UInt{1};
    constexpr UInt exponent_mask = (std::numeric_limits<UInt>::max() >> 1) & ~fraction_mask;
    // Preserve infinities and NaN payloads without floating-point arithmetic,
    // including signaling NaNs. Rounding finite values may overflow to infinity.
    if ((value & exponent_mask) == exponent_mask) return value;
    return transform_word(value, n, operation);
}

template <typename T>
bool copy_unmodified(const real_info_params<T> &params, const T *input,
                     size_t n_elem, int n, T *output) {
    validate_bit_count<T>(n, params.fixed_point);
    if (n_elem == 0) return true;
    if (input == nullptr || output == nullptr) {
        throw fixed_point_error("Bit operations require valid input and output arrays");
    }
    if (copy_constant_fixed_point(params, input, n_elem, output)) return true;
    if (n != 0) return false;
    if (params.fixed_point) {
        for (size_t i = 0; i < n_elem; ++i) {
            // Validate the input range, but avoid an encode/decode round trip.
            const FP<T> checked(params.fixed_range, input[i]);
            (void)checked;
        }
    }
    if (input != output) std::memmove(output, input, n_elem * sizeof(T));
    return true;
}

} // namespace real_info_detail

/**
 * @brief Unsigned native/fixed-point storage with shared low-bit operations.
 * Native counts are limited to the mantissa; fixed-point counts to 63 bits.
 */
template <typename T>
class bits {
    static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>);
    const real_info_params<T> &params_;

    void apply(int n, real_info_detail::bit_operation operation) {
        real_info_detail::validate_bit_count<T>(n, params_.fixed_point);
        if (n == 0 || (params_.fixed_point && params_.fixed_range.constant)) return;
        if (params_.fixed_point) {
            const auto maximum = static_cast<uint64_t>(FP<T>::max_code);
            if (bit_set_64 > maximum) throw fixed_point_error("Fixed-point code is outside the initialized mapping");
            // Setting, half-shaving, or rounding near the upper endpoint can
            // exceed Q. Saturate in unsigned arithmetic before decoding.
            bit_set_64 = std::min(real_info_detail::transform_word(bit_set_64, n, operation), maximum);
        } else if constexpr (std::is_same_v<T, float>) {
            bit_set_32 = real_info_detail::transform_native<T>(bit_set_32, n, operation);
        } else {
            bit_set_64 = real_info_detail::transform_native<T>(bit_set_64, n, operation);
        }
    }

public:
    uint32_t bit_set_32 = 0;
    uint64_t bit_set_64 = 0;

    bits(const real_info_params<T> &params, T value) : params_(params) {
        if (params_.fixed_point) {
            bit_set_64 = static_cast<uint64_t>(FP<T>(params_.fixed_range, value).value);
        } else if constexpr (std::is_same_v<T, float>) {
            bit_set_32 = std::bit_cast<uint32_t>(value);
        } else {
            bit_set_64 = std::bit_cast<uint64_t>(value);
        }
    }

    T convert_back() {
        if (params_.fixed_point) {
            FP<T> value(params_.fixed_range);
            value.value = std::bit_cast<int64_t>(bit_set_64);
            return value.convert_back();
        }
        if constexpr (std::is_same_v<T, float>) return std::bit_cast<T>(bit_set_32);
        else return std::bit_cast<T>(bit_set_64);
    }

    void shave(int n) { apply(n, real_info_detail::bit_operation::shave); }
    void set(int n) { apply(n, real_info_detail::bit_operation::set); }
    void halfshave(int n) { apply(n, real_info_detail::bit_operation::halfshave); }
    void round(int n) { apply(n, real_info_detail::bit_operation::round); }
};

/** 
 * @brief Shave the lowest `n` bits from each element in an array.
 * @param A Input array of type `T`.
 * @param n_elem Number of elements in the array.
 * @param n Number of least-significant bits to shave.
 * @param shaved Output array of shaved values.
 */
template <typename T>
void shave_bits(const real_info_params<T> &RI_params, T *A, size_t n_elem, int n, T *shaved) {
	if (real_info_detail::copy_unmodified(RI_params, A, n_elem, n, shaved)) return;
	for (size_t i = 0; i < n_elem; i++) {
		bits<T> b(RI_params, A[i]);
		uint64_t a;
		if (RI_params.fixed_point) {
			a = b.bit_set_64;
		}
		b.shave(n);
		if (RI_params.fixed_point) {
			FP<T> diff_fp(RI_params.fixed_range);
			diff_fp.value = static_cast<int64_t>(a - b.bit_set_64);
			const double diff = diff_fp.convert_back_no_offset();
			const double result = std::clamp(static_cast<double>(A[i]) - diff,
				static_cast<double>(RI_params.fixed_range.minimum),
				static_cast<double>(RI_params.fixed_range.maximum));
			shaved[i] = static_cast<T>(result);
		} else {
			shaved[i] = b.convert_back();
		}
	}
}

/** 
 * @brief Half-shave the lowest `n` bits from each element in an array.
 * @param RI_params Real information parameters.
 * @param A Input array of type `T`.
 * @param n_elem Number of elements in the array.
 * @param n Number of least-significant bits to half-shave.
 * @param halfshaved Output array of half-shaved values.
 */
template <typename T>
void halfshave_bits(const real_info_params<T> &RI_params, T *A, size_t n_elem, int n, T *halfshaved) {
	if (real_info_detail::copy_unmodified(RI_params, A, n_elem, n, halfshaved)) return;
	for (size_t i = 0; i < n_elem; i++) {
		bits<T> b(RI_params, A[i]);
		b.halfshave(n);
		halfshaved[i] = b.convert_back();
	}
}

/** 
 * @brief Bitwise round the `n` least-significant bits from each element in an array.
 * @param RI_params Real information parameters.
 * @param A Input array of type `T`.
 * @param n_elem Number of elements in the array.
 * @param n Number of least-significant bits to round.
 * @param bit_rounding_result Output array of rounded values.
 */
template <typename T>
void bit_rounding_bits(const real_info_params<T> &RI_params, T *A, size_t n_elem, int n, T *bit_rounding_result) {
	if (real_info_detail::copy_unmodified(RI_params, A, n_elem, n, bit_rounding_result)) return;
	for (size_t i = 0; i < n_elem; i++) {
		bits<T> b(RI_params, A[i]);
		b.round(n);
		bit_rounding_result[i] = b.convert_back();
	}
}

/** 
 * @brief Set the lowest `n` bits for each element in an array.
 * @param RI_params Real information parameters.
 * @param A Input array of type `T`.
 * @param n_elem Number of elements in the array.
 * @param n Number of least-significant bits to set.
 * @param set Output array of set values.
 */
template <typename T>
void set_bits(const real_info_params<T> &RI_params, T *A, size_t n_elem, int n, T *set) {
	if (real_info_detail::copy_unmodified(RI_params, A, n_elem, n, set)) return;
	for (size_t i = 0; i < n_elem; i++) {
		bits<T> b(RI_params, A[i]);
		b.set(n);
		set[i] = b.convert_back();
	}
}

/**
 * @brief Alternate setting and shaving the lowest `n` bits of array elements.
 * @param RI_params Real information parameters.
 * @param A Input array of type `T`.
 * @param n_elem Number of elements in the array.
 * @param n Number of least-significant bits to modify.
 * @param groomed Output array, with even-indexed elements set and odd ones shaved.
 */
template <typename T>
void groom_bits(const real_info_params<T> &RI_params, T *A, size_t n_elem, int n, T *groomed) {
	if (real_info_detail::copy_unmodified(RI_params, A, n_elem, n, groomed)) return;
	for (size_t i = 0; i < n_elem; i++) {
		if (i % 2 == 0) {
			set_bits<T>(RI_params, A + i, 1, n, groomed + i);
		} else {
			shave_bits<T>(RI_params, A + i, 1, n, groomed + i);
		}
	}
}

static int bit_pair_count_calls = 0;

static int n_exponent_bits_float = 8;
static int n_exponent_bits_double = 11;
static constexpr int n_joint_counts = 4;

// Per-bit work is bounded by the larger of the native and fixed-point widths,
// independent of the input length (including 63-bit fixed-point float inputs).
template <typename T, typename Value, size_t CountsPerBit = 1>
using bit_metric_buffer = std::array<Value,
    CountsPerBit * std::max(sizeof(T), sizeof(int64_t)) * CHAR_BIT>;

namespace real_info_detail {

template <typename T>
void require_analysis_input(const T *input, size_t n_elem) {
    if (input == nullptr || n_elem == 0) {
        throw fixed_point_error("Information analysis requires a nonempty array.");
    }
}

} // namespace real_info_detail

/**
 * @brief Counts the number of bits set to 1 in a specific bit position across an array using the `bits` class.
 * @param RI_params Real information parameters.
 * @param A Input array of type `T`.
 * @param n_elem Number of elements in the array.
 * @param i Bit position to count (0 for least significant bit).
 * @return The count of bits set to 1; zero if no samples remain after masking.
 */
template <typename T>
size_t bit_count_bits(const real_info_params<T> &RI_params, T *A, size_t n_elem, int i, size_t *non_zero) {
	size_t count = 0;
	*non_zero = 0;
	for (size_t j = 0; j < n_elem; j++) {
		if (RI_params.mask_zero && std::abs(A[j]) < EPS) {
			continue;
		}
		bits<T> b(RI_params, A[j]);
		if (RI_params.fixed_point) {
			uint64_t mask = 0x0000000000000001;
			mask <<= i;
			count += (b.bit_set_64 & mask) ? 1 : 0;
		} else if (std::is_same_v<T, float>) {
			uint32_t mask = 0x00000001;
			mask <<= i;
			count += (b.bit_set_32 & mask) ? 1 : 0;
		} else if (std::is_same_v<T, double>) {
			uint64_t mask = 0x0000000000000001;
			mask <<= i;
			count += (b.bit_set_64 & mask) ? 1 : 0;
		}
		(*non_zero)++;
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
void bit_count_bits(const real_info_params<T> &RI_params, T *A, size_t n_elem, size_t *c, size_t *non_zero) {
	int n_bits = sizeof(A[0]) * CHAR_BIT;
	int start_bit, end_bit;
	get_start_end_bits<T>(RI_params, n_bits, &start_bit, &end_bit);
	for (int i = start_bit; i < end_bit; i++) {
		c[i - start_bit] = bit_count_bits(RI_params, A, n_elem, i, non_zero);
	}
}

/**
 * @brief Computes the Shannon entropy of a probability distribution.
 * @param A Input array of probabilities (must sum to 1).
 * @param n_elem Number of elements in the array.
 * @return The Shannon entropy of the distribution represented by A.
 */
inline double entropy(double *A, size_t n_elem) {
	double H = 0;

	for (size_t i = 0; i < n_elem; i++) {
		if (A[i] > 0) {
			H += -A[i] * std::log2(A[i]);
		}
	}

	return H;
}

namespace real_info_detail {

struct paired_bit_information {
    double mutual = 0;
    double entropy_a = 0;
    double entropy_b = 0;
};

// Derive every paired metric from the same retained rows (00, 01, 10, 11).
// No retained rows means zero measured information, without dividing by zero.
inline paired_bit_information information_from_joint_counts(const size_t *counts, size_t n_samples) {
    paired_bit_information result;
    if (n_samples == 0) return result;

    const double n = static_cast<double>(n_samples);
    double pa[2] = {static_cast<double>(counts[0] + counts[1]) / n,
                    static_cast<double>(counts[2] + counts[3]) / n};
    double pb[2] = {static_cast<double>(counts[0] + counts[2]) / n,
                    static_cast<double>(counts[1] + counts[3]) / n};
    result.entropy_a = entropy(pa, 2);
    result.entropy_b = entropy(pb, 2);
    for (int a = 0; a < 2; ++a) {
        for (int b = 0; b < 2; ++b) {
            if (counts[2 * a + b] == 0) continue;
            const double p = static_cast<double>(counts[2 * a + b]) / n;
            result.mutual += p * std::log2(p / (pa[a] * pb[b]));
        }
    }
    // Roundoff in the sum must not turn an independent bit into negative information.
    result.mutual = std::max(0.0, result.mutual);
    return result;
}

} // namespace real_info_detail

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

  std::vector<T> A_sorted(A, A + n_elem);
  std::sort(A_sorted.begin(), A_sorted.end());
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
  return H;
}

/**
 * @brief Helper function to compute pair counts at a particular position in two numbers represented as bits.
 * @param RI_params Real information parameters.
 * @param a First number represented as bits.
 * @param b Second number represented as bits.
 * @param j Bit position to compare.
 * @param p Output array to store the counts of each bit pair.
 */
template <typename T>
void bit_pair_count_bits(const real_info_params<T> &RI_params, const bits<T> &a, const bits<T> &b, int j, int *p) {
	int bit_a, bit_b;

	if (RI_params.fixed_point) {
		uint64_t mask = 0x0000000000000001;
		mask <<= j;
		bit_a = (a.bit_set_64 & mask) ? 1 : 0;
		bit_b = (b.bit_set_64 & mask) ? 1 : 0;
	} else if (std::is_same_v<T, float>) {
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
 * @param RI_params Real information parameters.
 * @param A Pointer to the first array.	
 * @param B Pointer to the second array.
 * @param n_elem The number of elements in the arrays A and B.
 * @param pair_counts Pointer to an array to store the counts of each bit pair.
 * Counts are overwritten. With masking enabled, skip only pairs with both values masked.
 */
template <typename T>
void bit_pair_count_bits(const real_info_params<T> &RI_params, T *A, T *B, size_t n_elem, size_t *pair_counts, size_t *non_zero) {
	bit_pair_count_calls += 1;
	const int n_bits = sizeof(T) * CHAR_BIT;
	int start_bit, end_bit;
	get_start_end_bits<T>(RI_params, n_bits, &start_bit, &end_bit);
	int total_bits = end_bit - start_bit;
	std::fill_n(pair_counts, total_bits * n_joint_counts, size_t{0});
	*non_zero = 0;
	size_t n_non_zero = 0;
	#pragma omp parallel for reduction(+:pair_counts[:total_bits * n_joint_counts], n_non_zero)
	for (size_t i = 0; i < n_elem; i++) {
		if (RI_params.mask_zero && std::abs(A[i]) < EPS && std::abs(B[i]) < EPS) {
			continue;
		}
		bits<T> a(RI_params, A[i]);
		bits<T> b(RI_params, B[i]);
		int pair_count[4];
		for (int j = start_bit; j < end_bit; j++) {
			bit_pair_count_bits(RI_params, a, b, j, pair_count);
			pair_counts[(j - start_bit) * n_joint_counts] += pair_count[0];
			pair_counts[(j - start_bit) * n_joint_counts + 1] += pair_count[1];
			pair_counts[(j - start_bit) * n_joint_counts + 2] += pair_count[2];
			pair_counts[(j - start_bit) * n_joint_counts + 3] += pair_count[3];
		}
		++n_non_zero;
	}
	*non_zero = n_non_zero;
}

/**
 * @brief Computes the entropy of each bit position in an array.
 *
 * @tparam T The type of the elements in the input array.
 * @param A Pointer to the array of numeric values.
 * @param n_elem The number of elements in the array A.
 * @param H Pointer to an array to store the entropy of each bit position.
 * Nonempty arrays with every value masked have zero entropy at every bit.
 */
template <typename T>
void bit_count_entropy(const real_info_params<T> &RI_params, T *A, size_t n_elem, double *H) {
	real_info_detail::require_analysis_input(A, n_elem);
	const int n_bits = sizeof(A[0]) * CHAR_BIT;
	int start_bit, end_bit;
	get_start_end_bits<T>(RI_params, n_bits, &start_bit, &end_bit);
	bit_metric_buffer<T, size_t> c{};
	size_t n_non_zero;
	bit_count_bits(RI_params, A, n_elem, c.data(), &n_non_zero); // bits class implementation
	if (n_non_zero == 0) {
		std::fill_n(H, end_bit - start_bit, 0.0);
		return;
	}

	for (int i = start_bit; i < end_bit; i++) {
		double p = static_cast<double>(c[i - start_bit])/static_cast<double>(n_non_zero);
		double p_array[2] = {p, 1 - p};
		H[i - start_bit] = entropy(p_array, 2);
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
 * @param RI_params Real information parameters.
 * Nonempty arrays with no eligible pairs have zero mutual information at every bit.
 */
template <typename T>
void mutual_information(const real_info_params<T> &RI_params, T *A, T *B, size_t n_elem, double *info) {
	real_info_detail::require_analysis_input(A, n_elem);
	real_info_detail::require_analysis_input(B, n_elem);
	const int n_bits = sizeof(A[0]) * CHAR_BIT;
	int start_bit, end_bit;
	get_start_end_bits<T>(RI_params, n_bits, &start_bit, &end_bit);
	int total_bits = end_bit - start_bit;

	bit_metric_buffer<T, size_t, n_joint_counts> pair_counts{};
	size_t n_non_zero;
	bit_pair_count_bits(RI_params, A, B, n_elem, pair_counts.data(), &n_non_zero); // bits class implementation
	for (int i = 0; i < total_bits; i++) {
		info[i] = real_info_detail::information_from_joint_counts(
			pair_counts.data() + i * n_joint_counts, n_non_zero).mutual;
	}
}

// XXX: doing 1D for now, depends how multidimensional arrays are stored. 
// If using mdspan things would be straightforward, however there is poor support (for now).
// If using pointer, need to do some fancy(ish) indexing.
/** 
 * @brief Compute real information stored in given array by evaluating mutual information of array with a shifted version of itself
 * 
 * @tparam T Type of elements in array
 * @param RI_params Real information parameters
 * @param A Pointer to array
 * @param n_elem number of elements in A
 * @param info Pointer to array storing real information
 */
template <typename T>
void bitwise_real_information(const real_info_params<T> &RI_params, T *A, size_t n_elem, double *info) {
	real_info_detail::require_analysis_input(A, n_elem);
	if (n_elem == 1) {
		int start_bit, end_bit;
		get_start_end_bits<T>(RI_params, sizeof(T) * CHAR_BIT, &start_bit, &end_bit);
		std::fill_n(info, end_bit - start_bit, 0.0);
		return;
	}
	T *B = A + 1;
	mutual_information(RI_params, A, B, n_elem - 1, info);
}

/**
 * @brief Computes amount of redundant information between two arrays
 *
 * @tparam The types of arrays A and B
 * @param RI_params Real information parameters
 * @param A Pointer to first array 
 * @param B Pointer to first array 
 * @param n_elem Number of elements in each array 
 * @param R Pointer to array to store redundancy for each bit position 
 * No eligible pairs gives zero redundancy. Constant bits give one only if equal.
 */
template <typename T>
void redundancy(const real_info_params<T> &RI_params, T *A, T *B, size_t n_elem, double *R) {
	real_info_detail::require_analysis_input(A, n_elem);
	real_info_detail::require_analysis_input(B, n_elem);
	const int n_bits = sizeof(A[0]) * CHAR_BIT;
	int start_bit, end_bit;
	get_start_end_bits<T>(RI_params, n_bits, &start_bit, &end_bit);
	int total_bits = end_bit - start_bit;
	bit_metric_buffer<T, size_t, n_joint_counts> pair_counts{};
	size_t n_non_zero;
	bit_pair_count_bits(RI_params, A, B, n_elem, pair_counts.data(), &n_non_zero);

	for (int i = 0; i < total_bits; i++) {
		R[i] = 0.0;
		if (n_non_zero == 0) continue;
		const size_t *counts = pair_counts.data() + i * n_joint_counts;
		const auto metrics = real_info_detail::information_from_joint_counts(counts, n_non_zero);
		if (metrics.entropy_a + metrics.entropy_b > 0) {
			// The mathematical bounds are [0, 1]; clamp floating-point roundoff.
			R[i] = std::clamp(2 * metrics.mutual / (metrics.entropy_a + metrics.entropy_b), 0.0, 1.0);
		} else if (metrics.mutual == 0 && metrics.entropy_a == 0 && metrics.entropy_b == 0) {
			// Both bits are constant: equal constants are lossless, different ones are not.
			R[i] = counts[0] + counts[3] == n_non_zero ? 1.0 : 0.0;
		}
	}

}

/**
 *  @brief Computes real information preserved in array B with respect to A.
 *
 *  @tparam T Type of elements in arrays A and B
 *  @param RI_params Real information parameters
 *  @param A Pointer to original array
 *  @param B Pointer to shaved array
 *  @param n_elem Number of elements in each array
 *  @return Fraction of real information preserved in B with respect to A
 */
template <typename T>
double preserved_information_template(const real_info_params<T> &RI_params, T *A, T *B, size_t n_elem) {
	if (B == nullptr) {
		throw fixed_point_error("Information preservation requires a valid comparison array.");
	}
	if (is_constant_array(A, n_elem)) {
		// A constant reference has zero information: use an explicit lossless
		// convention instead of a ratio with a zero denominator.
		return std::memcmp(A, B, n_elem * sizeof(T)) == 0 ? 1.0 : 0.0;
	}
    real_info_params<T> float_params = RI_params;
	float_params.fixed_point = false;
	const int n_bits = sizeof(A[0]) * CHAR_BIT;
	int start_bit, end_bit;
	get_start_end_bits<T>(float_params, n_bits, &start_bit, &end_bit);
	int total_bits = end_bit - start_bit;
	bit_metric_buffer<T, double> R{};
	bit_metric_buffer<T, double> I{};

	redundancy(float_params, A, B, n_elem, R.data());
	bitwise_real_information(float_params, A, n_elem, I.data());

	double P = 0;
	double total_information = 0;
	for (int i = 0; i < total_bits; i++) {
		P += R[i] * I[i] ;
		total_information += I[i];
	}

	if (total_information == 0) {
		// Extend the constant-reference convention to any zero-information input.
		return std::memcmp(A, B, n_elem * sizeof(T)) == 0 ? 1.0 : 0.0;
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
 * @param RI_params Real information parameters
 */
template <typename T>
int pick_bits_to_shave_template(const real_info_params<T> &RI_params, T *A, size_t n_elem, double tolerance, int nbits_old) {
  if (is_constant_array(A, n_elem)) return 0;
  int max_bits = sizeof(T) * CHAR_BIT; 
  std::vector<T> scratch(n_elem);
  T *s = scratch.data();
  int start_bit;
  if (RI_params.fixed_point) {
	max_bits = sizeof(int64_t) * CHAR_BIT;
	start_bit = max_bits - 1;
  } else if constexpr (std::is_same_v<T, float>) { 
    start_bit = max_bits - n_exponent_bits_float - 1;
  } else if constexpr (std::is_same_v<T, double>) {
    start_bit = max_bits - n_exponent_bits_double - 1;
  }
  if (nbits_old == 0) {
    /* Start with shaving more bits rather than less */
    for (int i = start_bit; i >= 0; i--) {
		shave_bits<T>(RI_params, A, n_elem, i, s); // bits class implementation
		auto info = preserved_information_template<T>(RI_params, A, s, n_elem);
    	if (info > tolerance) {
    		return i;
    	}
    }
  } else {
    int nbits = nbits_old;
    	
    /* Check if we should retain more bits first */
	shave_bits<T>(RI_params, A, n_elem, nbits, s); // bits class implementation
	auto info = preserved_information_template<T>(RI_params, A, s, n_elem);
    if (info <= tolerance) {
      for (int i = nbits - 1; i >= 0; i--) {
		shave_bits<T>(RI_params, A, n_elem, i, s); // bits class implementation
		auto info = preserved_information_template<T>(RI_params, A, s, n_elem);
      	if (info > tolerance) {
      		return i;
      	}
      }
    } else {
  	/* Maybe we aren't throwing out enough. Increase number of bits thrown out */
  	for (int i = nbits + 1; i < max_bits; i++) {
		shave_bits<T>(RI_params, A, n_elem, i, s);
		auto info = preserved_information_template<T>(RI_params, A, s, n_elem);
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
 * @param RI_params Real information parameters.
 * @return The maximum number of least significant bits that can be cleared while
 *         preserving the exponent bits and sign bit.
 */
template <typename T>
int get_max_shave_bits(const real_info_params<T> &RI_params, int requested_cap = -1) {
  if (requested_cap != -1) {
    real_info_detail::validate_bit_count<T>(requested_cap, RI_params.fixed_point);
    return requested_cap;
  }
  int max_bits = sizeof(T) * CHAR_BIT; 
  if (RI_params.fixed_point) {
	max_bits = sizeof(int64_t) * CHAR_BIT;
	return max_bits - 1;
  } else if constexpr (std::is_same_v<T, float>) { 
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
void binary_search(const real_info_params<T> &RI_params, T *A, size_t n_elem, double tolerance, int &start_bit, int high_bit, int low_bit, T *shaved) {
    // Allocate once per search and reuse both candidates as the bounds narrow.
    std::vector<T> s_high(n_elem);
    std::vector<T> s_low(n_elem);
    for (;;) {
        shave_bits<T>(RI_params, A, n_elem, start_bit, s_high.data());
        shave_bits<T>(RI_params, A, n_elem, start_bit - 1, s_low.data());
        const auto info_high = preserved_information_template<T>(RI_params, A, s_high.data(), n_elem);
        const auto info_low = preserved_information_template<T>(RI_params, A, s_low.data(), n_elem);
        // Preserve the existing decisions and output behavior; boundary repairs
        // and the monotonicity requirement are tracked separately in step 8.
        if (info_high < tolerance && info_low >= tolerance) {
            std::copy(s_low.begin(), s_low.end(), shaved);
            --start_bit;
            return;
        } else if (info_low < tolerance) {
            high_bit = start_bit;
            start_bit = (start_bit - low_bit) / 2 + low_bit;
            if (start_bit == 0) {
                std::copy(s_low.begin(), s_low.end(), shaved);
                return;
            }
        } else if (info_high >= tolerance) {
            low_bit = start_bit;
            start_bit = (high_bit - start_bit) / 2 + start_bit;
            if (start_bit == low_bit) {
                std::copy(s_high.begin(), s_high.end(), shaved);
                return;
            }
        } else {
            std::cout << "binary search error start bit " << start_bit << " old bits shaved " << start_bit << " info high " << info_high << " info low " << info_low << " tolerance " << tolerance << "\n";
            return;
        }
    }
}

/**
 * @brief Picks the number of bits to shave from elements in array A such that the preserved information is above a given tolerance.
 *
 * @tparam T The type of the elements in the input array.
 * @param RI_params Real information parameters.
 * @param A Pointer to the input array.
 * @param n_elem The number of elements in the input array.
 * @param tolerance The minimum fraction of preserved information required.
 * @param nbits_old The previous number of bits shaved (0 if starting fresh).
 * @return The number of bits to shave from each element in the array.
 * @param config_file The path to the configuration file.
 */
template <typename T>
int pick_bits_to_shave_binary_search_template(const real_info_params<T> &RI_params, T *A, size_t n_elem, double tolerance, int nbits_old, char *config_file, T *shaved) {
	if (is_constant_array(A, n_elem)) {
		if (shaved == nullptr) {
			throw fixed_point_error("Bit selection requires a valid output array.");
		}
		if (A != shaved) std::copy_n(A, n_elem, shaved);
		return 0;
	}
	int max_shave_bits = get_max_shave_bits<T>(RI_params);

  // XXX: Check why do we need this?
  if (RI_params.mask_zero) {
    size_t n_non_zero = 0;
	for (size_t i = 0; i < n_elem; i++) {
	  if (std::abs(A[i]) >= EPS) {
		n_non_zero++;
	  }
	}
    if (n_non_zero == 0) {
      return 0;
    }
  }

  bit_pair_count_calls = 0;
  int start_bit;
  if (nbits_old == 0) {
	start_bit = max_shave_bits/2;
	binary_search<T>(RI_params, A, n_elem, tolerance, start_bit, max_shave_bits, 0, shaved);
  } else {
	start_bit = nbits_old;
	binary_search<T>(RI_params, A, n_elem, tolerance, start_bit, max_shave_bits, 0, shaved);
  }

  return start_bit;
}

/** Inclusive bounded search. A count and its candidate are committed together.
 * The information score is assumed monotone, as in the legacy binary selector.
 * Uncapped entry points above retain their existing selection behaviour.
 */
template <typename T>
int pick_bits_to_shave_binary_search_capped_template(
    const real_info_params<T> &params, T *input, size_t count, double tolerance,
    int previous, int requested_cap, T *output) {
    const int cap = get_max_shave_bits(params, requested_cap);
    if (count == 0 || input == nullptr || output == nullptr || input == output) {
        throw fixed_point_error("Capped selection requires nonempty input and a separate output array.");
    }
    if (!std::isfinite(tolerance) || tolerance < 0 || tolerance > 1) {
        throw fixed_point_error("Information tolerance must be finite and between 0 and 1.");
    }
    // Initialise every early-return path, using the representation's zero-bit
    // operation (which still performs the configured fixed-point conversion).
    shave_bits(params, input, count, 0, output);
    if (cap == 0 || is_constant_array(input, count)) return 0;
    if (params.mask_zero && std::none_of(input, input + count,
            [](T value) { return std::abs(value) >= EPS; })) return 0;

    std::vector<T> candidate(count);
    int low = 1, high = cap, selected = 0;
    int trial = std::clamp(previous > 0 ? previous : cap / 2, low, high);
    while (low <= high) {
        shave_bits(params, input, count, trial, candidate.data());
        const double preserved = preserved_information_template(params, input, candidate.data(), count);
        if (std::isfinite(preserved) && preserved >= tolerance) {
            selected = trial;
            std::copy(candidate.begin(), candidate.end(), output);
            low = trial + 1;
        } else {
            high = trial - 1;
        }
        if (low <= high) trial = low + (high - low) / 2;
    }
    return selected;
}

#endif
