#pragma once

#include <array>
#include <iomanip>
#include <optional>
#include <sstream>


namespace ruckig {

//! Result type of the OTGs update function
enum Result {
    Working = 0, ///< The trajectory is calculated normally
    Finished = 1, ///< Trajectory has reached its final position
    Error = -1, ///< Unclassified error
    ErrorInvalidInput = -100, ///< Error in the input parameter
    ErrorTrajectoryDuration = -101, ///< The trajectory duration exceed the numeral limits
    ErrorExecutionTimeCalculation = -110, ///< Error during the extremel time calculation (Step 1)
    ErrorSynchronizationCalculation = -111, ///< Error during the synchronization calculation (Step 2)
};


enum class Interface {
    Position, ///< Position-control: Full control over the entire kinematic state (Default)
    Velocity, ///< Velocity-control: Ignores the current position, target position, and velocity limits
};

enum class Synchronization {
    Time, ///< Always synchronize the DoFs to reach the target at the same time (Default)
    TimeIfNecessary, ///< Synchronize only when necessary (e.g. for non-zero target velocity or acceleration)
    None, ///< Calculate every DoF independently
};

enum class DurationDiscretization {
    Continuous, ///< Every trajectory duration is allowed (Default)
    Discrete, ///< The trajectory duration must be a multiple of the control cycle
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

    Interface interface {Interface::Position};
    Synchronization synchronization {Synchronization::Time};
    DurationDiscretization duration_discretization {DurationDiscretization::Continuous};

    Vector current_position, current_velocity {}, current_acceleration {};
    Vector target_position, target_velocity {}, target_acceleration {};
    Vector max_velocity, max_acceleration, max_jerk;
    std::optional<Vector> min_velocity, min_acceleration;

    std::array<bool, DOFs> enabled;
    std::optional<double> minimum_duration;

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
            || min_acceleration != rhs.min_acceleration
            || interface != rhs.interface
            || synchronization != rhs.synchronization
            || duration_discretization != rhs.duration_discretization
        );
    }

    std::string to_string() const {
        std::stringstream ss;
        ss << "\ninp.current_position = [" << this->join(current_position) << "]\n";
        ss << "inp.current_velocity = [" << this->join(current_velocity) << "]\n";
        ss << "inp.current_acceleration = [" << this->join(current_acceleration) << "]\n";
        ss << "inp.target_position = [" << this->join(target_position) << "]\n";
        ss << "inp.target_velocity = [" << this->join(target_velocity) << "]\n";
        ss << "inp.target_acceleration = [" << this->join(target_acceleration) << "]\n";
        ss << "inp.max_velocity = [" << this->join(max_velocity) << "]\n";
        ss << "inp.max_acceleration = [" << this->join(max_acceleration) << "]\n";
        ss << "inp.max_jerk = [" << this->join(max_jerk) << "]\n";
        if (min_velocity) {
            ss << "inp.min_velocity = [" << this->join(min_velocity.value()) << "]\n";
        }
        if (min_acceleration) {
            ss << "inp.min_acceleration = [" << this->join(min_acceleration.value()) << "]\n";
        }
        return ss.str();
    }
};

} // namespace ruckig
