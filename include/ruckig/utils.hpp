#pragma once

#include <string>
#include <sstream>
#include <iomanip>
#include <type_traits>

namespace ruckig {
 
    template<class Vector>
    std::string join(const Vector& array) {
        std::ostringstream ss;
        for (size_t i = 0; i < array.size(); ++i) {
            if (i) ss << ", ";
            ss << std::setprecision(16) << array[i];
        }
        return ss.str();
    }

} // namespace ruckig
