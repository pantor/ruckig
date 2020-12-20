#pragma once


namespace ruckig {

inline double Power(double v, int e) {
    return std::pow(v, e);
}

inline double Power(double v, double e) {
    return std::pow(v, e);
}

inline double Sqrt(double v) {
    return std::sqrt(v);
}

inline double Abs(double v) {
    return std::abs(v);
}

} // namespace ruckig
