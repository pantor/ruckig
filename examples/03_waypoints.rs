use ruckig::{InputParameter, OutputParameter, Result, Ruckig};

fn main() {
    let control_cycle = 0.01;
    let dofs = 3;
    let max_number_of_waypoints = 10;  // for memory allocation

    // Create instances: Ruckig as well as input and output parameters
    let mut ruckig = Ruckig::new(dofs, control_cycle, max_number_of_waypoints);
    let mut input = InputParameter::new(dofs);
    let mut output = OutputParameter::new(dofs, max_number_of_waypoints);

    // Set input parameters
    input.current_position = vec![0.2, 0.0, -0.3];
    input.current_velocity = vec![0.0, 0.2, 0.0];
    input.current_acceleration = vec![0.0, 0.6, 0.0];

    input.intermediate_positions = vec![
        vec![1.4, -1.6, 1.0],
        vec![-0.6, -0.5, 0.4],
        vec![-0.4, -0.35, 0.0],
        vec![0.8, 1.8, -0.1]
    ];

    input.target_position = vec![0.5, 1.0, 0.0];
    input.target_velocity = vec![0.2, 0.0, 0.3];
    input.target_acceleration = vec![0.0, 0.1, -0.1];

    input.max_velocity = vec![1.0, 2.0, 1.0];
    input.max_acceleration = vec![3.0, 2.0, 2.0];
    input.max_jerk = vec![6.0, 10.0, 20.0];

    // Generate the trajectory within the control loop
    println!("t | position");
    while ruckig.update(&mut input, &mut output) == Result::Working {
        let pos: Vec<String> = output.new_position.iter().map(|x| format!("{x:.4}")).collect();
        println!("{:.4} | [{}]", output.time, pos.join(", "));

        output.pass_to_input(&mut input);
    }

    println!("Trajectory duration: {} [s].", output.trajectory().get_duration());
}

