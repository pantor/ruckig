#include <ruckig/ruckig.hpp>


namespace ruckig {

inline double v_at_t(double v0, double a0, double j, double t) {
    return v0 + t * (a0 + j * t / 2);
}

void Brake::acceleration_brake(double v0, double a0, double vMax, double aMax, double jMax, std::array<double, 2>& t_brake, std::array<double, 2>& j_brake) {
    j_brake[0] = -jMax;

    double t_to_a_max = (a0 - aMax) / jMax;
    double v_at_a_max = v_at_t(v0, a0, -jMax, t_to_a_max);

    t_brake[0] = t_to_a_max + eps;

    double t_to_a_zero = a0 / jMax;
    double v_at_a_zero = v_at_t(v0, a0, -jMax, t_to_a_zero);

    if ((v_at_a_zero > vMax && jMax > 0) || (v_at_a_zero < vMax && jMax < 0)) {
        t_brake[0] = t_to_a_zero + eps;
        v_at_a_max = v_at_a_zero;
    }

    if ((v_at_a_max < -vMax && jMax > 0) || (v_at_a_max > -vMax && jMax < 0)) {
        double t_to_v_max_while_a_max = -(v_at_a_max + vMax)/aMax;
        double t_to_v_max_in_reverse_j_direction = -aMax/(2*jMax) + (vMax - v_at_a_max)/aMax;
        t_brake[1] = std::min(t_to_v_max_while_a_max, t_to_v_max_in_reverse_j_direction);
        
    } else if ((v_at_a_max > vMax && jMax > 0) || (v_at_a_max < vMax && jMax < 0)) {
        double t_to_other_a_max = (a0 + aMax) / jMax - eps;
        double t_to_v_max = a0/jMax + std::sqrt(a0*a0 + 2 * jMax * (v0 - vMax)) / std::abs(jMax);
        double t_to_v_max_brake_for_other = a0/jMax + std::sqrt(a0*a0/2 + jMax * (v0 + vMax)) / std::abs(jMax);

        if (t_to_v_max_brake_for_other < t_to_other_a_max && t_to_v_max_brake_for_other < t_to_v_max) {
            t_brake[0] = t_to_v_max_brake_for_other - eps;

        } else if (t_to_other_a_max < t_to_v_max) {
            double v_at_a_other_a_max = v_at_t(v0, a0, -jMax, t_to_other_a_max);
            double t_to_v_max_while_a_max = (v_at_a_other_a_max - vMax)/aMax;
            double t_to_v_max_brake_for_other_a_max = -aMax/(2*jMax) + (vMax + v_at_a_other_a_max)/aMax;

            t_brake[0] = t_to_other_a_max - eps;
            t_brake[1] = std::min(t_to_v_max_while_a_max, t_to_v_max_brake_for_other_a_max);
        
        } else {
            t_brake[0] = t_to_v_max + eps;
        }
    }
}

void Brake::velocity_brake(double v0, double a0, double vMax, double aMax, double jMax, std::array<double, 2>& t_brake, std::array<double, 2>& j_brake) {
    j_brake[0] = -jMax;
    double t_to_a_max = (a0 + aMax)/jMax;
    double t_to_v_max_in_j_direction = a0/jMax + std::sqrt(a0*a0 + 2 * jMax * (v0 - vMax)) / std::abs(jMax);
    double t_to_v_max_in_reverse_j_direction = a0/jMax + std::sqrt(a0*a0/2 + jMax * (v0 + vMax)) / std::abs(jMax);

    if (t_to_a_max < t_to_v_max_in_j_direction && t_to_a_max < t_to_v_max_in_reverse_j_direction) {
        double v_at_a_max = v_at_t(v0, a0, -jMax, t_to_a_max);
        double t_to_v_max_while_a_max = (v_at_a_max - vMax)/aMax;
        double t_to_v_max_in_reverse_j_direction = -aMax/(2*jMax) + (v_at_a_max + vMax)/aMax;

        t_brake[0] = t_to_a_max;
        t_brake[1] = std::min(t_to_v_max_while_a_max, t_to_v_max_in_reverse_j_direction);
    } else {
        t_brake[0] = std::min(t_to_v_max_in_j_direction, t_to_v_max_in_reverse_j_direction);
    }

    t_brake[0] = std::max(t_brake[0] - eps, 0.0);
}

void Brake::get_brake_trajectory(double v0, double a0, double vMax, double vMin, double aMax, double jMax, std::array<double, 2>& t_brake, std::array<double, 2>& j_brake) {
    t_brake[0] = 0.0;
    t_brake[1] = 0.0;
    j_brake[0] = 0.0;
    j_brake[1] = 0.0;

    const double a0_sq_jMax = std::pow(a0, 2)/(2 * jMax) - vMax;

    if (a0 > aMax) {
        acceleration_brake(v0, a0, vMax, aMax, jMax, t_brake, j_brake);

    } else if (a0 < -aMax) {
        acceleration_brake(v0, a0, vMin, -aMax, -jMax, t_brake, j_brake);

    } else if ((v0 > vMax && (a0 > 0 || a0_sq_jMax < v0)) || (a0 > 0 && a0_sq_jMax > -v0)) {
        velocity_brake(v0, a0, vMax, aMax, jMax, t_brake, j_brake);

    } else if ((v0 < vMin && (a0 < 0 || a0_sq_jMax > v0)) || (a0 < 0 && a0_sq_jMax > v0)) {
        velocity_brake(v0, a0, vMin, -aMax, -jMax, t_brake, j_brake);
    }
}

} // namespace ruckig
