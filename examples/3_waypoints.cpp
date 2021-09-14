#include <iostream>

#include <ruckig/ruckig.hpp>


using namespace ruckig;

int main() {
    // Create instances: the Ruckig OTG as well as input and output parameters
    Ruckig<3> otg {0.01};  // control cycle
    InputParameter<3> input;
    OutputParameter<3> output;

    // Set input parameters
    input.current_position = {0.0, 0.0, 0.0};
    input.current_velocity = {0.0, 0.0, 0.0};
    input.current_acceleration = {0.0, 0.0, 0.0};

    input.intermediate_positions = {
        {1.4, -1.6, 1.0},
        {-0.6, -0.5, 0.4},
        {-0.4, -0.35, 0.0},
        {0.8, 1.8, -0.1}
    };

    input.target_position = {1.0, 1.0, 1.0};
    input.target_velocity = {0.0, 0.0, 0.0};
    input.target_acceleration = {0.0, 0.0, 0.0};

    input.max_velocity = {1.0, 1.0, 1.0};
    input.max_acceleration = {2.0, 2.0, 2.0};
    input.max_jerk = {10.0, 10.0, 10.0};

    std::cout << "t | p1 | p2 | p3" << std::endl;
    while (otg.update(input, output) == Result::Working) {
        auto& p = output.new_position;
        std::cout << output.time << " " << p[0] << " " << p[1] << " " << p[2] << std::endl;

        input.current_position = output.new_position;
        input.current_velocity = output.new_velocity;
        input.current_acceleration = output.new_acceleration;
    }

    std::cout << "Reached target position in " << output.trajectory.get_duration() << " [s]." << std::endl;
    std::cout << "Calculation in " << output.calculation_duration << " [µs]." << std::endl;
}
