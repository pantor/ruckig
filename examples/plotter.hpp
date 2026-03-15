#pragma once

#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>


//! Join a vector for pretty printing (e.g. to std::cout)
template<class Vector>
inline std::string pretty_print(const Vector& array) {
    std::ostringstream ss;
    for (size_t i = 0; i < array.size(); ++i) {
        if (i) ss << ", ";
        ss << array[i];
    }
    return ss.str();
}
