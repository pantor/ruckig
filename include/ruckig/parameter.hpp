#pragma once

#include <array>
#include <iomanip>
#include <optional>
#include <sstream>


namespace ruckig {

//! Result type of the OTGs update function
enum Result {
    Working = 0,
    Finished = 1,
    Error = -1,
    ErrorInvalidInput = -100,
    ErrorTrajectoryDuration = -101,
    ErrorExecutionTimeCalculation = -110,
    ErrorSynchronizationCalculation = -111,
};


//! Input type of the OTG
template<size_t DOFs>
class InputParameter {
    template<class T>
    static std::string join(const T& array) {
        std::ostringstream ss;
        for (size_t i = 0; i < DOFs; ++i) {
            if (i) ss << ", ";
            ss << std::setprecision(15) << array[i];
        }
        return ss.str();
    }

public:
    using Vector = std::array<double, DOFs>;
    static constexpr size_t degrees_of_freedom {DOFs};

    enum class Type {
        Position,
        Velocity,
    } type {Type::Position};

    Vector current_position, current_velocity {}, current_acceleration {};
    Vector target_position, target_velocity {}, target_acceleration {};
    Vector max_velocity, max_acceleration, max_jerk;

    std::array<bool, DOFs> enabled;
    std::optional<double> minimum_duration;
    std::optional<Vector> min_velocity;

    InputParameter() {
        std::fill(enabled.begin(), enabled.end(), true);
    }

    bool operator!=(const InputParameter<DOFs>& rhs) const {
        return (
            current_position != rhs.current_position
            || current_velocity != rhs.current_velocity
            || current_acceleration != rhs.current_acceleration
            || target_position != rhs.target_position
            || target_velocity != rhs.target_velocity
            || target_acceleration != rhs.target_acceleration
            || max_velocity != rhs.max_velocity
            || max_acceleration != rhs.max_acceleration
            || max_jerk != rhs.max_jerk
            || enabled != rhs.enabled
            || minimum_duration != rhs.minimum_duration
            || min_velocity != rhs.min_velocity
            || type != rhs.type
        );
    }

    std::string to_string() const {
        std::stringstream ss;
        ss << "\ninp.current_position = [" << join(current_position) << "]\n";
        ss << "inp.current_velocity = [" << join(current_velocity) << "]\n";
        ss << "inp.current_acceleration = [" << join(current_acceleration) << "]\n";
        ss << "inp.target_position = [" << join(target_position) << "]\n";
        ss << "inp.target_velocity = [" << join(target_velocity) << "]\n";
        ss << "inp.target_acceleration = [" << join(target_acceleration) << "]\n";
        ss << "inp.max_velocity = [" << join(max_velocity) << "]\n";
        ss << "inp.max_acceleration = [" << join(max_acceleration) << "]\n";
        ss << "inp.max_jerk = [" << join(max_jerk) << "]\n";
        if (min_velocity) {
            ss << "inp.min_velocity = [" << join(min_velocity.value()) << "]\n";
        }
        return ss.str();
    }
};


//! Output type of the OTG
template<size_t DOFs>
struct OutputParameter {
    using Vector = std::array<double, DOFs>;
    static constexpr size_t degrees_of_freedom {DOFs};

    Vector new_position, new_velocity, new_acceleration;

    double duration; // [s]
    bool new_calculation {false};
    double calculation_duration; // [µs]

    Vector independent_min_durations; // [s]
};

} // namespace ruckig
