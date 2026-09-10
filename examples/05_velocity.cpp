#include <ruckig/ruckig.hpp>

#include "plotter.hpp"


using namespace ruckig;

int main() {
    const size_t degrees_of_freedom = 3;
    const double control_cycle = 0.01;

    // Create instances: the Ruckig trajectory generator as well as input and output parameters
    Ruckig<DynamicDOFs> ruckig(degrees_of_freedom, control_cycle);
    InputParameter<DynamicDOFs> input(degrees_of_freedom);
    OutputParameter<DynamicDOFs> output(degrees_of_freedom);

    // Set input parameters and velocity control interface
    input.control_interface = ControlInterface::Velocity;

    input.current_position = {0.0, 0.0, 0.5};
    input.current_velocity = {3.0, -2.2, -0.5};
    input.current_acceleration = {0.0, 2.5, -0.5};

    input.target_velocity = {0.0, -0.5, -1.5};
    input.target_acceleration = {0.0, 0.0, 0.5};

    input.max_acceleration = {3.0, 2.0, 1.0};
    input.max_jerk = {6.0, 6.0, 4.0};

    // Generate the trajectory within the control loop
    std::cout << "t | position" << std::endl;
    while (ruckig.update(input, output) == Result::Working) {
        std::cout << output.time << " | " << pretty_print(output.new_position) << std::endl;

        output.pass_to_input(input);
    }

    std::cout << "Trajectory duration: " << output.trajectory.get_duration() << " [s]." << std::endl;
}
