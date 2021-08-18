#pragma once

#include <array>
#include <iomanip>
#include <optional>
#include <sstream>
#include <type_traits>
#include <vector>


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


enum class ControlInterface {
    Position, ///< Position-control: Full control over the entire kinematic state (Default)
    Velocity, ///< Velocity-control: Ignores the current position, target position, and velocity limits
};

enum class Synchronization {
    Phase, ///< Phase synchronize the DoFs when possible, else fallback to "Time" strategy
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
    template<class T> using Vector = typename std::conditional<DOFs >= 1, std::array<T, DOFs>, std::vector<T>>::type;

    static std::string join(const Vector<double>& array) {
        std::ostringstream ss;
        for (size_t i = 0; i < array.size(); ++i) {
            if (i) ss << ", ";
            ss << std::setprecision(16) << array[i];
        }
        return ss.str();
    }

    void initialize() {
        std::fill(current_velocity.begin(), current_velocity.end(), 0.0);
        std::fill(current_acceleration.begin(), current_acceleration.end(), 0.0);
        std::fill(target_velocity.begin(), target_velocity.end(), 0.0);
        std::fill(target_acceleration.begin(), target_acceleration.end(), 0.0);
        std::fill(enabled.begin(), enabled.end(), true);
    }

public:
    size_t degrees_of_freedom;

    ControlInterface control_interface {ControlInterface::Position};
    Synchronization synchronization {Synchronization::Time};
    DurationDiscretization duration_discretization {DurationDiscretization::Continuous};

    //! Current state
    Vector<double> current_position, current_velocity, current_acceleration;

    //! Target state
    Vector<double> target_position, target_velocity, target_acceleration;

    //! Kinematic constraints
    Vector<double> max_velocity, max_acceleration, max_jerk;
    std::optional<Vector<double>> min_velocity, min_acceleration;

    //! Intermediate waypoints (not yet used)
    std::vector<Vector<double>> intermediate_positions;

    //! Is the DoF considered for calculation?
    Vector<bool> enabled;

    //! Per-DoF control_interface (overwrites global synchronization, not yet used)
    // std::optional<Vector<ControlInterface>> per_dof_control_interface;

    //! Per-DoF synchronization (overwrites global synchronization, not yet used)
    // std::optional<Vector<Synchronization>> per_dof_synchronization;

    //! Optional minimum trajectory duration
    std::optional<double> minimum_duration;

    //! Optional duration [s] after which the trajectory calculation is (softly) interrupted (not yet used)
    std::optional<double> interrupt_calculation_duration;

    template <size_t D = DOFs, typename std::enable_if<D >= 1, int>::type = 0>
    InputParameter(): degrees_of_freedom(DOFs) {
        initialize();
    }

    template <size_t D = DOFs, typename std::enable_if<D == 0, int>::type = 0>
    InputParameter(size_t dofs): degrees_of_freedom(dofs) {
        current_position.resize(dofs);
        current_velocity.resize(dofs);
        current_acceleration.resize(dofs);
        target_position.resize(dofs);
        target_velocity.resize(dofs);
        target_acceleration.resize(dofs);
        max_velocity.resize(dofs);
        max_acceleration.resize(dofs);
        max_jerk.resize(dofs);
        enabled.resize(dofs);

        initialize();
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
            || intermediate_positions != rhs.intermediate_positions
            || enabled != rhs.enabled
            || minimum_duration != rhs.minimum_duration
            || min_velocity != rhs.min_velocity
            || min_acceleration != rhs.min_acceleration
            || control_interface != rhs.control_interface
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
        if (!intermediate_positions.empty()) {
            ss << "inp.intermediate_positions = [\n";
            for (auto p: intermediate_positions) {
                ss << "    [" << this->join(p) << "],\n";
            }
            ss << "]\n";
        }
        return ss.str();
    }
};

} // namespace ruckig
