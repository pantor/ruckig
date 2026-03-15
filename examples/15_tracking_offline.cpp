// Only with Ruckig Pro

#include <cmath>

#include <ruckig/trackig.hpp>

#include "plotter.hpp"


using namespace ruckig;

// Create the target state signal
TargetState<1> model_half_sinus(double t, double ramp_vel=1.2) {
    TargetState<1> target;
    if (t < 2.5) { // [s]
        target.position[0] = std::sin(ramp_vel * t);
        target.velocity[0] = ramp_vel * std::cos(ramp_vel * t);
        target.acceleration[0] = -ramp_vel * ramp_vel * std::sin(ramp_vel * t);
    } else {
        target.position[0] = 0.0;
        target.velocity[0] = 0.0;
        target.acceleration[0] = 0.0;
    }
    return target;
}


int main() {
    // Create instances: the Trackig trajectory generator as well as input and output parameters
    Trackig<1> trackig(0.01);  // control cycle
    InputParameter<1> input;

    // Set input parameters
    input.current_position = {0.0};
    input.current_velocity = {0.0};
    input.current_acceleration = {0.0};

    input.max_velocity = {4.0};
    input.max_acceleration = {5.0};
    input.max_jerk = {15.0};

    // trackig.look_ahead_cycles = 32;  // look ahead in the future, the transition between the sinus and the plateau should get smoother

    // Pre-generate the trajectory
    std::vector< TargetState<1> > target_states;
    for (size_t t = 0; t < 400; ++t) {
        target_states.push_back(model_half_sinus(trackig.delta_time * t));
    }

    // Calculate the full trajectory offline
    std::vector< OutputParameter<1> > output_states = trackig.calculate_trajectory(target_states, input);

    // Generate the trajectory following the target state
    std::cout << "target | follow" << std::endl;
    for (size_t t = 0; t < std::min<size_t>(target_states.size(), output_states.size()); t += 1) {
        std::cout << pretty_print(target_states[t].position) << " | " << pretty_print(output_states[t].new_position) << std::endl;
    }
}
