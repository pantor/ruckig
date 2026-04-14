#pragma once

#include <memory>

#include <ruckig/ruckig.hpp>


using Ruckig  = ruckig::Ruckig<0>;
using InputParameter = ruckig::InputParameter<0>;
using OutputParameter = ruckig::OutputParameter<0>;
using Trajectory = ruckig::Trajectory<0>;


static inline std::vector<double> to_vec(const double* data, size_t len) { return {data, data + len}; }


// Trajectory
std::unique_ptr<Trajectory> traj_new_direct(size_t dofs) { return std::make_unique<Trajectory>(dofs); }
std::unique_ptr<Trajectory> traj_new(size_t dofs, size_t waypoints) { return std::make_unique<Trajectory>(dofs, waypoints); }
size_t traj_get_dofs(const Trajectory& traj) { return traj.degrees_of_freedom; }
double traj_get_duration(const Trajectory& traj) { return traj.get_duration(); }
void traj_at_time(const Trajectory& traj, double time, double* pos, double* vel, double* acc, size_t len) {
    std::vector<double> p(len), v(len), a(len);
    traj.at_time(time, p, v, a);
    std::copy(p.begin(), p.end(), pos);
    std::copy(v.begin(), v.end(), vel);
    std::copy(a.begin(), a.end(), acc);
}

// InputParameter
std::unique_ptr<InputParameter> input_new(size_t dofs) { return std::make_unique<InputParameter>(dofs); };
int32_t input_get_control_interface(const InputParameter& input) { return static_cast<int32_t>(input.control_interface); }
void input_set_control_interface(InputParameter& input, int32_t v) { input.control_interface = static_cast<ruckig::ControlInterface>(v); }
int32_t input_get_synchronization(const InputParameter& input) { return static_cast<int32_t>(input.synchronization); }
void input_set_synchronization(InputParameter& input, int32_t v) { input.synchronization = static_cast<ruckig::Synchronization>(v); }
void input_set_current_position(InputParameter& input, const double* data, size_t len) { input.current_position = to_vec(data, len); }
const std::vector<double>& input_get_current_position(const InputParameter& input) { return input.current_position; }
void input_set_current_velocity(InputParameter& input, const double* data, size_t len) { input.current_velocity = to_vec(data, len); }
const std::vector<double>& input_get_current_velocity(const InputParameter& input) { return input.current_velocity; }
void input_set_current_acceleration(InputParameter& input, const double* data, size_t len) { input.current_acceleration = to_vec(data, len); }
const std::vector<double>& input_get_current_acceleration(const InputParameter& input) { return input.current_acceleration; }
void input_set_target_position(InputParameter& input, const double* data, size_t len) { input.target_position = to_vec(data, len); }
const std::vector<double>& input_get_target_position(const InputParameter& input) { return input.target_position; }
void input_set_target_velocity(InputParameter& input, const double* data, size_t len) { input.target_velocity = to_vec(data, len); }
const std::vector<double>& input_get_target_velocity(const InputParameter& input) { return input.target_velocity; }
void input_set_target_acceleration(InputParameter& input, const double* data, size_t len) { input.target_acceleration = to_vec(data, len); }
const std::vector<double>& input_get_target_acceleration(const InputParameter& input) { return input.target_acceleration; }
void input_set_max_velocity(InputParameter& input, const double* data, size_t len) { input.max_velocity = to_vec(data, len); }
const std::vector<double>& input_get_max_velocity(const InputParameter& input) { return input.max_velocity; }
void input_set_max_acceleration(InputParameter& input, const double* data, size_t len) { input.max_acceleration = to_vec(data, len); }
const std::vector<double>& input_get_max_acceleration(const InputParameter& input) { return input.max_acceleration; }
void input_set_max_jerk(InputParameter& input, const double* data, size_t len) { input.max_jerk = to_vec(data, len); }
const std::vector<double>& input_get_max_jerk(const InputParameter& input) { return input.max_jerk; }
void input_set_max_position(InputParameter& input, const double* data, size_t len) { input.max_position = to_vec(data, len); }
const std::vector<double>& input_get_max_position(const InputParameter& input) { return input.max_position; }
void input_set_min_position(InputParameter& input, const double* data, size_t len) { input.min_position = to_vec(data, len); }
const std::vector<double>& input_get_min_position(const InputParameter& input) { return input.min_position; }
void input_set_intermediate_positions(InputParameter& input, const double* data, size_t dofs, size_t count) {
    input.intermediate_positions.resize(count);
    for (size_t r = 0; r < count; ++r) {
        input.intermediate_positions[r] = to_vec(data + r * dofs, dofs);
    }
}
void input_clear_intermediate_positions(InputParameter& input) { input.intermediate_positions.clear(); }
size_t input_get_intermediate_positions(const InputParameter& input, double* buf, size_t capacity) {
    size_t written = 0;
    for (auto& row : input.intermediate_positions) {
        if (written + row.size() > capacity) break;
        std::copy(row.begin(), row.end(), buf + written);
        written += row.size();
    }
    return written;
}

// OutputParameter
std::unique_ptr<OutputParameter> output_new_direct(size_t dofs) { return std::make_unique<OutputParameter>(dofs); };
std::unique_ptr<OutputParameter> output_new(size_t dofs, size_t waypoints) { return std::make_unique<OutputParameter>(dofs, waypoints); };
const Trajectory& output_get_trajectory(const OutputParameter& output) { return output.trajectory; }
const std::vector<double>& output_get_new_position(const OutputParameter& output) { return output.new_position; }
const std::vector<double>& output_get_new_velocity(const OutputParameter& output) { return output.new_velocity; }
const std::vector<double>& output_get_new_acceleration(const OutputParameter& output) { return output.new_acceleration; }
double output_get_time(const OutputParameter& output) { return output.time; };
size_t output_get_new_section(const OutputParameter& output) { return output.new_section; }
bool output_get_new_calculation(const OutputParameter& output) { return output.new_calculation; }
double output_get_calculation_duration(const OutputParameter& output) { return output.calculation_duration; }
void output_pass_to_input(const OutputParameter& output, InputParameter& input) { output.pass_to_input(input); }

// Ruckig
std::unique_ptr<Ruckig> ruckig_new_direct_and_offline(size_t dofs) { return std::make_unique<Ruckig>(dofs); }
std::unique_ptr<Ruckig> ruckig_new_direct(size_t dofs, double delta_time) { return std::make_unique<Ruckig>(dofs, delta_time); }
std::unique_ptr<Ruckig> ruckig_new_offline(size_t dofs, size_t waypoints) { return std::make_unique<Ruckig>(dofs, waypoints); }
std::unique_ptr<Ruckig> ruckig_new(size_t dofs, double delta_time, size_t waypoints) { return std::make_unique<Ruckig>(dofs, delta_time, waypoints); }
int ruckig_update(Ruckig& ruckig, const InputParameter& input, OutputParameter& output) { return static_cast<int>(ruckig.update(input, output)); }
int ruckig_calculate(Ruckig& ruckig, const InputParameter& input, Trajectory& trajectory) { return static_cast<int>(ruckig.calculate(input, trajectory)); }
void ruckig_reset(Ruckig& ruckig) { ruckig.reset(); }
