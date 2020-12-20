#include <ruckig/ruckig.hpp>


namespace ruckig {

inline double v_at_t(double v0, double a0, double j, double t) {
    return v0 + a0 * t + j * std::pow(t, 2) / 2;
}

void RuckigStep1::get_brake_trajectory(double v0, double a0, double vMax, double aMax, double jMax, std::array<double, 2>& t_brake, std::array<double, 2>& j_brake) {
    t_brake[0] = 0.0;
    t_brake[1] = 0.0;
    j_brake[0] = 0.0;
    j_brake[1] = 0.0;

    const double a0_sq_jMax = std::pow(a0, 2)/(2 * jMax) - vMax;

    if (a0 > aMax) {
        j_brake[0] = -jMax;

        double t_to_a_max = (a0 - aMax) / jMax;
        double v_at_a_max = v_at_t(v0, a0, -jMax, t_to_a_max);

        t_brake[0] = t_to_a_max + 2e-15;

        double t_to_a_zero = a0 / jMax;
        double v_at_a_zero = v_at_t(v0, a0, -jMax, t_to_a_zero);

        if (v_at_a_zero > vMax) {
            t_brake[0] = t_to_a_zero + 2e-15;
            v_at_a_max = v_at_a_zero;
        }

        if (v_at_a_max < -vMax) {
            double t_to_v_max_while_a_max = -(v_at_a_max + vMax)/aMax;
            double t_to_v_max_in_reverse_j_direction = (-std::pow(aMax, 2) - 2 * jMax * v_at_a_max + 2 * jMax * vMax)/(2*aMax*jMax);
            t_brake[1] = std::min(t_to_v_max_while_a_max, t_to_v_max_in_reverse_j_direction);
            
        } else if (v_at_a_max > vMax) {
            double t_to_other_a_max = (a0 + aMax) / jMax - 2e-15;
            double t_to_v_max = -(-a0 - std::sqrt(std::pow(a0,2) + 2 * jMax * v0 - 2 * jMax * vMax))/jMax;
            double t_to_v_max_brake_for_other = (2 * a0 + std::sqrt(2) * std::sqrt(std::pow(a0, 2) + 2 * jMax * (v0 + vMax)))/(2 * jMax);

            // std::cout << t_to_other_a_max << " " << t_to_v_max << " " << t_to_v_max_brake_for_other << std::endl;

            if (t_to_v_max_brake_for_other < t_to_other_a_max && t_to_v_max_brake_for_other < t_to_v_max) {
                t_brake[0] = t_to_v_max_brake_for_other - 2e-15;

            } else if (t_to_other_a_max < t_to_v_max) {
                double v_at_a_other_a_max = v_at_t(v0, a0, -jMax, t_to_other_a_max);
                double t_to_v_max_while_a_max = (v_at_a_other_a_max - vMax)/aMax;
                double t_to_v_max_brake_for_other_a_max = -(std::pow(aMax, 2) - 2 * jMax * (v_at_a_max + vMax))/(2 * aMax * jMax);

                t_brake[0] = t_to_other_a_max;
                t_brake[1] = std::min(t_to_v_max_while_a_max, t_to_v_max_brake_for_other_a_max);
            
            } else {
                t_brake[0] = t_to_v_max + 2e-15;
            }
        }

    } else if (a0 < -aMax) {
        j_brake[0] = jMax;

        double t_to_a_max = -(a0 + aMax) / jMax;
        double v_at_a_max = v_at_t(v0, a0, jMax, t_to_a_max);

        t_brake[0] = t_to_a_max + 2e-15;

        double t_to_a_zero = -a0 / jMax;
        double v_at_a_zero = v_at_t(v0, a0, jMax, t_to_a_zero);
        if (v_at_a_zero < -vMax) {
            t_brake[0] = t_to_a_zero + 2e-15;
            v_at_a_max = v_at_a_zero;
        }

        if (v_at_a_max > vMax) {
            double t_to_v_max_while_a_max = (v_at_a_max - vMax)/aMax;
            double t_to_v_max_in_reverse_j_direction = (-std::pow(aMax, 2) + 2 * jMax * v_at_a_max + 2 * jMax * vMax)/(2*aMax*jMax);
            t_brake[1] = std::min(t_to_v_max_while_a_max, t_to_v_max_in_reverse_j_direction);
            // {a1 == -aMax, v1 == v0 - aMax * t1, -vMax == v1 + a1 * t2 + jMax / 2 * t2^2, 0 == a1 + jMax t2}, {t1, t2, a1, v1}

        } else if (v_at_a_max < -vMax) {
            double t_to_other_a_max = -(a0 - aMax) / jMax - 2e-15;
            double t_to_v_max = (-a0 + std::sqrt(std::pow(a0,2) - 2 * jMax * (v0 + vMax)))/jMax;

            if (t_to_other_a_max < t_to_v_max) {
                double v_at_a_other_a_max = v_at_t(v0, a0, jMax, t_to_other_a_max);
                double t_to_v_max_while_a_max = -(v_at_a_other_a_max + vMax)/aMax;
                double t_to_v_max_brake_for_other_a_max = -(std::pow(aMax, 2) + 2 * jMax * (v_at_a_max - vMax))/(2 * aMax * jMax);
                
                t_brake[0] = t_to_other_a_max;
                t_brake[1] = std::min(t_to_v_max_while_a_max, t_to_v_max_brake_for_other_a_max);
            } else {
                t_brake[0] = t_to_v_max + 2e-15;
            }
        }

    } else if ((v0 > vMax && (a0 > 0 || a0_sq_jMax < v0)) || (a0 > 0 && a0_sq_jMax > -v0)) {
        j_brake[0] = -jMax;
        double t_to_a_max = (a0 + aMax)/jMax;
        double t_to_v_max_in_j_direction = (a0 + std::sqrt(std::pow(a0, 2) - 2 * jMax * (vMax - v0)))/jMax;
        double t_to_v_max_in_reverse_j_direction = (2 * a0 + std::sqrt(2*(std::pow(a0, 2) + 2 * jMax * (v0 + vMax))))/(2 * jMax);

        if (t_to_a_max < t_to_v_max_in_j_direction && t_to_a_max < t_to_v_max_in_reverse_j_direction) {
            double v_at_a_max = v_at_t(v0, a0, -jMax, t_to_a_max);
            double t_to_v_max_while_a_max = (v_at_a_max - vMax)/aMax;
            double t_to_v_max_in_reverse_j_direction = -(std::pow(aMax, 2) - 2 * jMax * (v_at_a_max + vMax))/(2 * aMax * jMax);

            t_brake[0] = t_to_a_max;
            t_brake[1] = std::min(t_to_v_max_while_a_max, t_to_v_max_in_reverse_j_direction);
        } else {
            t_brake[0] = std::min(t_to_v_max_in_j_direction, t_to_v_max_in_reverse_j_direction);
        }

        t_brake[0] = std::max(t_brake[0] - 2e-15, 0.0);

    } else if ((v0 < -vMax && (a0 < 0 || a0_sq_jMax > v0)) || (a0 < 0 && a0_sq_jMax > v0)) {
        j_brake[0] = jMax;
        double t_to_a_max = -(a0 - aMax)/jMax;
        double t_to_v_max_in_j_direction = (-a0 + std::sqrt(std::pow(a0, 2) + 2 * jMax * (-vMax - v0)))/jMax;
        double t_to_v_max_in_reverse_j_direction = (-2 * a0 + std::sqrt(2*(std::pow(a0, 2) - 2 * jMax * (v0 - vMax))))/(2 * jMax);

        if (t_to_a_max < t_to_v_max_in_j_direction && t_to_a_max < t_to_v_max_in_reverse_j_direction) {
            double v_at_a_max = v_at_t(v0, a0, jMax, t_to_a_max);
            double t_to_v_max_while_a_max = -(v_at_a_max + vMax)/aMax;
            double t_to_v_max_in_reverse_j_direction = -(std::pow(aMax, 2) + 2 * jMax * (v_at_a_max - vMax))/(2 * aMax * jMax);

            t_brake[0] = t_to_a_max;
            t_brake[1] = std::min(t_to_v_max_while_a_max, t_to_v_max_in_reverse_j_direction);
        } else {
            t_brake[0] = std::min(t_to_v_max_in_j_direction, t_to_v_max_in_reverse_j_direction);
        }

        t_brake[0] = std::max(t_brake[0] - 2e-15, 0.0);
    }

    // std::cout << t_brake[0] << " " << t_brake[1] << std::endl;
}

} // namespace ruckig
