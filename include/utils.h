#pragma once
#include <iostream>
#include <string>
#include <vector>

template<typename T>
void debug_cout(std::string field_name, const T &value) {
	// return;
	std::cout << "\"+++" << field_name << "\": " << value << "," << std::endl;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& arr) {
    os << "[";
	for (int i = 0; i < arr.size(); i++) {
		os << arr[i];
		if (i + 1 < arr.size()) {
			os << ", ";
		}
	}
	os << "]";

    return os;
}

