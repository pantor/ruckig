use ruckig::{InputParameter, OutputParameter, Result, Ruckig};

fn main() {
    // Create instances: Ruckig as well as input and output parameters
    let mut ruckig = Ruckig::new_direct(3, 0.01); // DoFs, control cycle
    let mut input = InputParameter::new(3);
    let mut output = OutputParameter::new_direct(3);

    // Set input parameters and velocity control interface
    input.control_interface = ruckig::ControlInterface::Velocity;

    input.current_position = vec![0.0, 0.0, 0.5];
    input.current_velocity = vec![3.0, -2.2, -0.5];
    input.current_acceleration = vec![0.0, 2.5, -0.5];

    input.target_velocity = vec![0.0, -0.5, -1.5];
    input.target_acceleration = vec![0.0, 0.0, 0.5];

    input.max_acceleration = vec![3.0, 2.0, 1.0];
    input.max_jerk = vec![6.0, 6.0, 4.0];

    // Generate the trajectory within the control loop
    println!("t | position");
    while ruckig.update(&mut input, &mut output) == Result::Working {
        let pos: Vec<String> = output.new_position.iter().map(|x| format!("{x:.4}")).collect();
        println!("{:.4} | [{}]", output.time, pos.join(", "));

        output.pass_to_input(&mut input);
    }

    println!("Trajectory duration: {} [s].", output.trajectory().get_duration());
}

