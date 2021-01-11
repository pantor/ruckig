#include <iostream>

#include <ruckig/parameter.hpp>
#include <ruckig/ruckig.hpp>


using namespace ruckig;

int main() {
    Ruckig<6> otg {0.001};
    InputParameter<6> input;
    OutputParameter<6> output;

    input.max_velocity = {1.2, 1.2, 1.2, 0.6, 0.6, 0.6};
    input.max_acceleration = {4.0, 4.0, 4.0, 1.5, 1.5, 1.5};
    input.max_jerk = {10.0, 10.0, 10.0, 4.0, 4.0, 4.0};

    input.current_position = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    input.current_velocity = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    input.target_position = {1.0, 0.5, 0.5, 0.0, -0.1, 0.2};

    std::cout << "t | p1 | p2 | p3 | p4 | p5 | p6" << std::endl;
    while (otg.update(input, output) == Result::Working) {
        auto& new_p = output.new_position;
        std::cout << otg.get_time() << " " << new_p[0] << " " << new_p[1] << " " << new_p[2] << " " << new_p[3] << " " << new_p[4] << " " << new_p[5] << std::endl;

        input.current_position = output.new_position;
        input.current_velocity = output.new_velocity;
        input.current_acceleration = output.new_acceleration;
    }

    std::cout << "Reached target position in " << output.duration << " [s]." << std::endl;
}
