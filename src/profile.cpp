#include <iomanip>

#include <ruckig/ruckig.hpp>


namespace ruckig {

std::string Profile::to_string() const {
    std::string result;
    switch (direction) {
        case Direction::UP: result += "UP_"; break;
        case Direction::DOWN: result += "DOWN_"; break;
    }
    switch (limits) {
        case Limits::ACC0_ACC1_VEL: result += "ACC0_ACC1_VEL"; break;
        case Limits::VEL: result += "VEL"; break;
        case Limits::ACC0: result += "ACC0"; break;
        case Limits::ACC1: result += "ACC1"; break;
        case Limits::ACC0_ACC1: result += "ACC0_ACC1"; break;
        case Limits::ACC0_VEL: result += "ACC0_VEL"; break;
        case Limits::ACC1_VEL: result += "ACC1_VEL"; break;
        case Limits::NONE: result += "NONE"; break;
    }
    return result; 
}

std::tuple<double, double, double> Profile::integrate(double t, double p0, double v0, double a0, double j) {
    const double p_new = p0 + t * (v0 + t * (a0 / 2 + t * j / 6));
    const double v_new = v0 + t * (a0 + t * j / 2);
    const double a_new = a0 + t * j;
    return {p_new, v_new, a_new};
}

} // namespace ruckig
