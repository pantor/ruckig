# Only with Ruckig Pro

from math import sin, cos

from ruckig import Trackig, TargetState, InputParameter


# Create the target state signal
def model_half_sinus(t, ramp_vel=1.2):
    target = TargetState(1)
    if t < 2.5:  # [s]
        target.position = [sin(ramp_vel * t)]
        target.velocity = [ramp_vel * cos(ramp_vel * t)]
        target.acceleration = [-ramp_vel * ramp_vel * sin(ramp_vel * t)]
    else:
        target.position = [0.0]
        target.velocity = [0.0]
        target.acceleration = [0.0]
    return target


if __name__ == '__main__':
    # Create instances: the Trackig OTG as well as input and output parameters
    inp = InputParameter(1)
    trackig = Trackig(inp.degrees_of_freedom, 0.01)

    # Set input parameters
    inp.current_position = [0.0]
    inp.current_velocity = [0.0]
    inp.current_acceleration = [0.0]

    inp.max_velocity = [4.0]
    inp.max_acceleration = [5.0]
    inp.max_jerk = [15.0]

    # trackig.look_ahead_cycles = 32  # look ahead in the future, the transition between the sinus and the plateau should get smoother

    # Pre-generate the trajectory
    target_states = [model_half_sinus(trackig.delta_time * t) for t in range(400)]

    # Calculate the full trajectory offline
    output_states = trackig.calculate_trajectory(target_states, inp)

    print('target | follow')
    for target_state, out in zip(target_states, output_states):
        print(
            '\t'.join([f'{p:0.3f}' for p in target_state.position] + [f'{p:0.3f}' for p in out.new_position]),
            f'in {out.calculation_duration:0.2f} [µs]',
        )

    steps = [trackig.delta_time * i for i in range(len(output_states))]
    target_list = [[ts.position[0], ts.velocity[0], ts.acceleration[0]] for ts in target_states]
    follow_list = [[out.new_position[0], out.new_velocity[0], out.new_acceleration[0]] for out in output_states]

    # Plot the trajectory
    # from pathlib import Path
    # project_path = Path(__file__).parent.parent.absolute()

    # import numpy as np
    # import matplotlib.pyplot as plt

    # follow_list = np.array(follow_list)
    # target_list = np.array(target_list)

    # target_list = np.pad(target_list, [(0, 1), (0, 0)], 'constant')  # Pad to match lengths

    # plt.ylabel(f'DoF 1')
    # plt.plot(steps, follow_list[:, 0], label='Follow Position')
    # plt.plot(steps, follow_list[:, 1], label='Follow Velocity', linestyle='dotted')
    # plt.plot(steps, follow_list[:, 2], label='Follow Acceleration', linestyle='dotted')
    # plt.plot(steps, target_list[:, 0], color='r', label='Target Position')
    # plt.grid(True)
    # plt.legend()

    # plt.savefig(project_path / 'examples' / '15_trajectory.pdf')
