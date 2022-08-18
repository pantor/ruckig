// Only with Ruckig Pro

#include <iostream>

#include <ruckig/trackig.hpp>


using namespace ruckig;

TargetState<1> model_ramp(size_t t, double ramp_vel=0.5, double ramp_pos=1.0) {
    TargetState<1> target;
    const bool on_ramp = t < ramp_pos / (0.01 * std::abs(ramp_vel));
    target.position[0] = on_ramp ? t * ramp_vel * 0.01 : ramp_pos;
    target.velocity[0] = on_ramp ? ramp_vel : 0.0;
    target.acceleration[0] = 0.0;
    return target;
}

int main() {
    // Create instances: the Ruckig OTG as well as input and output parameters
    Trackig<1> otg {0.01};  // control cycle
    InputParameter<1> input;
    OutputParameter<1> output;

    // Set input parameters
    input.current_position = {0.0};
    input.current_velocity = {0.0};
    input.current_acceleration = {0.0};

    input.max_velocity = {0.8};
    input.max_acceleration = {2.0};
    input.max_jerk = {5.0};

    // Generate the trajectory within the control loop
    std::cout << "target | follow" << std::endl;
    for (size_t t = 0; t < 500; t += 1) {
        auto target_state = model_ramp(t);
        auto res = otg.update(target_state, input, output);
        std::cout << target_state.position[0] << " " << output.new_position[0] << std::endl;

        output.pass_to_input(input);
    }
}
