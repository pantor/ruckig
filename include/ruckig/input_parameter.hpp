#pragma once

#include <array>
#include <iomanip>
#include <optional>
#include <sstream>
#include <type_traits>
#include <vector>

#include <ruckig/result.hpp>
#include <ruckig/utils.hpp>

namespace ruckig {

enum class ControlInterface {
    Position, ///< Position-control: Full control over the entire kinematic state (Default)
    Velocity, ///< Velocity-control: Ignores the current position, target position, and velocity limits
};

enum class Synchronization {
    Time, ///< Always synchronize the DoFs to reach the target at the same time (Default)
    TimeIfNecessary, ///< Synchronize only when necessary (e.g. for non-zero target velocity or acceleration)
    Phase, ///< Phase synchronize the DoFs when this is possible, else fallback to "Time" strategy. Phase synchronization will result a straight-line trajectory
    // PhaseOnly, ///< Always phase synchronize the DoFs (even when this is not time-optimal), else returns "ErrorNoPhaseSynchronization". Ruckig will then guarantee a straight-line trajectory
    None, ///< Calculate every DoF independently
};

enum class DurationDiscretization {
    Continuous, ///< Every trajectory synchronization duration is allowed (Default)
    Discrete, ///< The trajectory synchronization duration must be a multiple of the control cycle
};


//! Input type of Ruckig
template<size_t DOFs, template<class, size_t> class CustomVector = StandardVector>
class InputParameter {
    template<class T> using Vector = CustomVector<T, DOFs>;

    void initialize() {
        for (size_t dof = 0; dof < degrees_of_freedom; ++dof) {
            current_velocity[dof] = 0.0;
            current_acceleration[dof] = 0.0;
            target_velocity[dof] = 0.0;
            target_acceleration[dof] = 0.0;
            enabled[dof] = true;
        }
    }

    void resize(size_t dofs) {
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
    }

#if defined WITH_ONLINE_CLIENT
    void reserve(size_t max_number_of_waypoints) {
        intermediate_positions.reserve(max_number_of_waypoints);
    }
#endif

public:
    size_t degrees_of_freedom;

    ControlInterface control_interface {ControlInterface::Position};
    Synchronization synchronization {Synchronization::Time};
    DurationDiscretization duration_discretization {DurationDiscretization::Continuous};

    // Current state
    Vector<double> current_position, current_velocity, current_acceleration;

    // Target state
    Vector<double> target_position, target_velocity, target_acceleration;

    // Kinematic constraints
    Vector<double> max_velocity, max_acceleration, max_jerk;
    std::optional<Vector<double>> min_velocity, min_acceleration;

    //! Intermediate waypoints (only in Ruckig Pro)
    std::vector<Vector<double>> intermediate_positions;

    // Kinematic constraints for intermediate sections (between waypoints) (only in Ruckig Pro)
    std::optional<std::vector<Vector<double>>> per_section_max_velocity, per_section_max_acceleration, per_section_max_jerk;
    std::optional<std::vector<Vector<double>>> per_section_min_velocity, per_section_min_acceleration;

    // Positional constraints (only in Ruckig Pro)
    std::optional<Vector<double>> max_position, min_position;

    //! Is the DoF considered for calculation?
    Vector<bool> enabled;

    //! Per-DoF control_interface (overwrites global control_interface)
    std::optional<Vector<ControlInterface>> per_dof_control_interface;

    //! Per-DoF synchronization (overwrites global synchronization)
    std::optional<Vector<Synchronization>> per_dof_synchronization;

    //! Optional minimum trajectory duration
    std::optional<double> minimum_duration;

    //! Optional minimum trajectory duration for each intermediate sections (only in Ruckig Pro)
    std::optional<std::vector<double>> per_section_minimum_duration;

    //! Optional duration [µs] after which the trajectory calculation is (softly) interrupted (only in Ruckig Pro)
    //!
    //! The total calculation consists of a first iterative phase and a second fixed phase. The interrupt signal
    //! is applied to the iterative phase only, and the real-time capable (constant) second phase is computed
    //! afterwards. Therefore, the total calculation duration might exceed this interrupt signal by a constant offset,
    //! which should be considered (subtracted) here.
    std::optional<double> interrupt_calculation_duration;

    template <size_t D = DOFs, typename std::enable_if<(D >= 1), int>::type = 0>
    InputParameter(): degrees_of_freedom(DOFs) {
        initialize();
    }

    template <size_t D = DOFs, typename std::enable_if<(D == 0), int>::type = 0>
    InputParameter(size_t dofs): degrees_of_freedom(dofs) {
        resize(dofs);
        initialize();
    }

#if defined WITH_ONLINE_CLIENT
    template <size_t D = DOFs, typename std::enable_if<(D >= 1), int>::type = 0>
    InputParameter(size_t max_number_of_waypoints): degrees_of_freedom(DOFs) {
        reserve(max_number_of_waypoints);
        initialize();
    }

    template <size_t D = DOFs, typename std::enable_if<(D == 0), int>::type = 0>
    InputParameter(size_t dofs, size_t max_number_of_waypoints): degrees_of_freedom(dofs) {
        reserve(max_number_of_waypoints);
        resize(dofs);
        initialize();
    }
#endif

    bool operator!=(const InputParameter<DOFs, CustomVector>& rhs) const {
        return !(
            current_position == rhs.current_position
            && current_velocity == rhs.current_velocity
            && current_acceleration == rhs.current_acceleration
            && target_position == rhs.target_position
            && target_velocity == rhs.target_velocity
            && target_acceleration == rhs.target_acceleration
            && max_velocity == rhs.max_velocity
            && max_acceleration == rhs.max_acceleration
            && max_jerk == rhs.max_jerk
            && intermediate_positions == rhs.intermediate_positions
            && per_section_max_velocity == rhs.per_section_max_velocity
            && per_section_max_acceleration == rhs.per_section_max_acceleration
            && per_section_max_jerk == rhs.per_section_max_jerk
            && per_section_min_velocity == rhs.per_section_min_velocity
            && per_section_min_acceleration == rhs.per_section_min_acceleration
            && max_position == rhs.max_position
            && min_position == rhs.min_position
            && enabled == rhs.enabled
            && minimum_duration == rhs.minimum_duration
            && per_section_minimum_duration == rhs.per_section_minimum_duration
            && min_velocity == rhs.min_velocity
            && min_acceleration == rhs.min_acceleration
            && control_interface == rhs.control_interface
            && synchronization == rhs.synchronization
            && duration_discretization == rhs.duration_discretization
            && per_dof_control_interface == rhs.per_dof_control_interface
            && per_dof_synchronization == rhs.per_dof_synchronization
        );
    }

    std::string to_string() const {
        std::stringstream ss;
        ss << "\ninp.current_position = [" << join(current_position, degrees_of_freedom) << "]\n";
        ss << "inp.current_velocity = [" << join(current_velocity, degrees_of_freedom) << "]\n";
        ss << "inp.current_acceleration = [" << join(current_acceleration, degrees_of_freedom) << "]\n";
        ss << "inp.target_position = [" << join(target_position, degrees_of_freedom) << "]\n";
        ss << "inp.target_velocity = [" << join(target_velocity, degrees_of_freedom) << "]\n";
        ss << "inp.target_acceleration = [" << join(target_acceleration, degrees_of_freedom) << "]\n";
        ss << "inp.max_velocity = [" << join(max_velocity, degrees_of_freedom) << "]\n";
        ss << "inp.max_acceleration = [" << join(max_acceleration, degrees_of_freedom) << "]\n";
        ss << "inp.max_jerk = [" << join(max_jerk, degrees_of_freedom) << "]\n";
        if (min_velocity) {
            ss << "inp.min_velocity = [" << join(min_velocity.value(), degrees_of_freedom) << "]\n";
        }
        if (min_acceleration) {
            ss << "inp.min_acceleration = [" << join(min_acceleration.value(), degrees_of_freedom) << "]\n";
        }
        if (!intermediate_positions.empty()) {
            ss << "inp.intermediate_positions = [\n";
            for (auto p: intermediate_positions) {
                ss << "    [" << join(p, degrees_of_freedom) << "],\n";
            }
            ss << "]\n";
        }
        return ss.str();
    }
};

} // namespace ruckig
