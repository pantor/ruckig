#include <iostream>

#include <ruckig/ruckig.hpp>


using namespace ruckig;

int main() {
    // Create instances: the ruckig otg as well as input and output parameters
    Ruckig<6> otg {0.001};
    InputParameter<6> input;
    OutputParameter<6> output;

    // Set input parameters
    input.max_velocity = {1.2, 1.2, 1.2, 0.6, 0.6, 0.6};
    input.max_acceleration = {4.0, 4.0, 4.0, 1.5, 1.5, 1.5};
    input.max_jerk = {10.0, 10.0, 10.0, 4.0, 4.0, 4.0};

    input.current_position = {0.0, -0.1, 0.12, 0.0, 0.3, 0.05};
    input.current_velocity = {0.0, 0.0, 0.2, 0.0, 0.0, 0.0};
    input.current_acceleration = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    input.target_position = {1.0, 0.5, 0.5, 0.0, -0.1, 0.2};
    input.target_velocity = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    input.target_acceleration = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // Generate trajectory
    std::cout << "t | p1 | p2 | p3 | p4 | p5 | p6" << std::endl;
    while (otg.update(input, output) == Result::Working) {
        auto& p = output.new_position;
        std::cout << output.time << " " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << " " << p[4] << " " << p[5] << std::endl;

        input.current_position = output.new_position;
        input.current_velocity = output.new_velocity;
        input.current_acceleration = output.new_acceleration;
    }

    std::cout << "Reached target position in " << output.trajectory.get_duration() << " [s]." << std::endl;
}
