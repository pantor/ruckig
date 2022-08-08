#pragma once

#include <array>
#include <string>
#include <sstream>
#include <iomanip>
#include <tuple>
#include <type_traits>
#include <vector>

namespace ruckig {


//! Constant for indicating a dynamic (run-time settable) number of DoFs
constexpr static size_t DynamicDOFs {0};


//! Vector data type based on the C++ standard library
template<class T, size_t DOFs> using StandardVector = typename std::conditional<DOFs >= 1, std::array<T, DOFs>, std::vector<T>>::type;
template<class T, size_t DOFs, size_t SIZE> using StandardSizeVector = typename std::conditional<DOFs >= 1, std::array<T, SIZE>, std::vector<T>>::type;


//! Vector data type based on the Eigen matrix type. Eigen needs to be included seperately
#ifdef EIGEN_VERSION_AT_LEAST
#if EIGEN_VERSION_AT_LEAST(3,4,0)
    template<class T, size_t DOFs> using EigenVector = typename std::conditional<DOFs >= 1, Eigen::Vector<T, DOFs>, Eigen::Vector<T, Eigen::Dynamic>>::type;
#endif
#endif


template<class Vector>
inline std::string join(const Vector& array, size_t size) {
    std::ostringstream ss;
    for (size_t i = 0; i < size; ++i) {
        if (i) ss << ", ";
        ss << std::setprecision(16) << array[i];
    }
    return ss.str();
}


//! Integrate with constant jerk for duration t. Returns new position, new velocity, and new acceleration.
inline std::tuple<double, double, double> integrate(double t, double p0, double v0, double a0, double j) {
    return std::make_tuple(
        p0 + t * (v0 + t * (a0 / 2 + t * j / 6)),
        v0 + t * (a0 + t * j / 2),
        a0 + t * j
    );
}

} // namespace ruckig
