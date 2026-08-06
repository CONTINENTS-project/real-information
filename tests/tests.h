#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <type_traits>
#include <cctype>

const size_t n_elements = 1000;

void test_assert(bool condition, const char* file, int line, std::string message) {
	if (!condition) {
		std::cerr << "Assertion failed: " << message << " at " << file << ":" << line << "\n";
		std::cerr << "FAIL " << file << "\n";
		exit(EXIT_FAILURE);
	}
}

bool float_equals(float a, float b, float epsilon = 1e-6f) {
	return std::fabs(a - b) < epsilon;
}

bool double_equals(double a, double b, double epsilon = 1e-12f) {
	return std::fabs(a - b) < epsilon;
}

void read_input_data(const char *filename, float *data, size_t n_elem) {
	std::ifstream infile(filename);
	for (size_t i = 0; i < n_elem; i++) {
		infile >> data[i];
	}
}

std::string trim(const std::string &s) {
	size_t start = s.find_first_not_of(" \t\n\r");
	if (start == std::string::npos) {
		return "";
	}
	size_t end = s.find_last_not_of(" \t\n\r");
	return s.substr(start, end - start + 1);
}

size_t leading_ws_count(const std::string &s) {
	size_t i = 0;
	while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) {
		++i;
	}
	return i;
}

std::string strip_quotes(const std::string &s) {
	if (s.size() >= 2) {
		if ((s.front() == '"' && s.back() == '"') ||
		    (s.front() == '\'' && s.back() == '\'')) {
			return s.substr(1, s.size() - 2);
		}
	}
	return s;
}

std::vector<std::string> parse_inline_yaml_list(const std::string &list_text) {
	std::string s = trim(list_text);
	if (s.size() < 2 || s.front() != '[' || s.back() != ']') {
		throw std::runtime_error("Expected inline YAML list: [..]");
	}

	std::string inner = s.substr(1, s.size() - 2);
	std::vector<std::string> out;
	std::string token;
	bool in_single = false;
	bool in_double = false;

	for (char c : inner) {
		if (c == '\'' && !in_double) {
			in_single = !in_single;
			token.push_back(c);
			continue;
		}
		if (c == '"' && !in_single) {
			in_double = !in_double;
			token.push_back(c);
			continue;
		}

		if (c == ',' && !in_single && !in_double) {
			std::string t = strip_quotes(trim(token));
			out.push_back(t);
			token.clear();
			continue;
		}

		token.push_back(c);
	}

	if (!token.empty() || inner.find(',') != std::string::npos) {
		std::string t = strip_quotes(trim(token));
		out.push_back(t);
	}

	return out;
}

std::string read_yaml_value(const std::string &filename, const std::string &key) {
	std::ifstream infile(filename);
	if (!infile.is_open()) {
		throw std::runtime_error("Failed to open YAML file: " + filename);
	}

	std::string line;
	while (std::getline(infile, line)) {
		std::string no_comment = line.substr(0, line.find('#'));
		std::string content = trim(no_comment);
		if (content.empty()) {
			continue;
		}

		size_t colon_pos = content.find(':');
		if (colon_pos == std::string::npos) {
			continue;
		}

		std::string current_key = trim(content.substr(0, colon_pos));
		if (current_key != key) {
			continue;
		}

		std::string value = trim(content.substr(colon_pos + 1));
		if (value.size() >= 2) {
			if ((value.front() == '"' && value.back() == '"') ||
			    (value.front() == '\'' && value.back() == '\'')) {
				value = value.substr(1, value.size() - 2);
			}
		}

		return value;
	}

	throw std::runtime_error("Key not found in YAML file: " + key);
}

std::vector<std::string> read_yaml_list(const std::string &filename, const std::string &key) {
	std::ifstream infile(filename);
	if (!infile.is_open()) {
		throw std::runtime_error("Failed to open YAML file: " + filename);
	}

	std::string line;
	bool in_block_list = false;
	bool found_key = false;
	size_t key_indent = 0;
	std::vector<std::string> values;

	while (std::getline(infile, line)) {
		std::string no_comment = line.substr(0, line.find('#'));
		std::string content = trim(no_comment);
		if (content.empty()) {
			continue;
		}

		size_t indent = leading_ws_count(no_comment);

		if (!in_block_list) {
			size_t colon_pos = content.find(':');
			if (colon_pos == std::string::npos) {
				continue;
			}

			std::string current_key = trim(content.substr(0, colon_pos));
			if (current_key != key) {
				continue;
			}

			found_key = true;
			key_indent = indent;

			std::string rhs = trim(content.substr(colon_pos + 1));
			if (!rhs.empty()) {
				return parse_inline_yaml_list(rhs);
			}

			in_block_list = true;
			continue;
		}

		if (indent <= key_indent) {
			break;
		}

		std::string t = trim(no_comment);
		if (!t.empty() && t.front() == '-') {
			std::string item = strip_quotes(trim(t.substr(1)));
			values.push_back(item);
		} else {
			break;
		}
	}

	if (!found_key) {
		throw std::runtime_error("Key not found in YAML file: " + key);
	}

	return values;
}

template <typename T>
T yaml_convert_scalar(const std::string &value_str) {
	if constexpr (std::is_same_v<T, float>) {
		return std::stof(value_str);
	} else if constexpr (std::is_same_v<T, double>) {
		return std::stod(value_str);
	} else if constexpr (std::is_same_v<T, int>) {
		return std::stoi(value_str);
	} else if constexpr (std::is_same_v<T, size_t>) {
		return static_cast<size_t>(std::stoull(value_str));
	} else if constexpr (std::is_same_v<T, std::string>) {
		return value_str;
	} else {
		throw std::runtime_error("Unexpected type for YAML value conversion");
	}
}

template <typename T>
void read_yaml_value(const std::string &filename, const std::string &key, T *value) {
	std::string value_str = read_yaml_value(filename, key);
	*value = yaml_convert_scalar<T>(value_str);
}

template <typename T>
void read_yaml_list(const std::string &filename, const std::string &key, T *out_values, size_t n_elem) {
	std::vector<std::string> raw = read_yaml_list(filename, key);
	if (raw.size() != n_elem) {
		throw std::runtime_error(
			"YAML list size mismatch for key '" + key + "': expected "
			+ std::to_string(n_elem) + ", got " + std::to_string(raw.size())
		);
	}

	for (size_t i = 0; i < n_elem; ++i) {
		out_values[i] = yaml_convert_scalar<T>(raw[i]);
	}
}