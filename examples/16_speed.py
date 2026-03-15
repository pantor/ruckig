# Only with Ruckig Pro

from copy import copy

from ruckig import InputParameter, OutputParameter, Result, Ruckig


if __name__ == '__main__':
    # Create instances: the Ruckig OTG as well as input and output parameters
    otg = Ruckig(3, 0.01)  # DoFs, control cycle
    inp = InputParameter(3)
    out = OutputParameter(3)

    # Set input parameters
    inp.current_position = [0.0, 0.0, 0.5]
    inp.current_velocity = [0.0, -2.2, -0.5]
    inp.current_acceleration = [0.0, 2.5, -0.5]

    inp.target_position = [5.0, -2.0, -3.5]
    inp.target_velocity = [0.0, -0.5, -2.0]
    inp.target_acceleration = [0.0, 0.0, 0.5]

    inp.max_velocity = [3.0, 1.0, 3.0]
    inp.max_acceleration = [3.0, 2.0, 1.0]
    inp.max_jerk = [4.0, 3.0, 2.0]

    # We want to break to a paused state and re-accelerate to the normal trajectory
    # using Ruckig's speed control feature.
    # In this example, we have the four phases: 1. start, 2. brake, 3. accelerate, 4. end.
    phase = 'start'
    speed_change_duration = 1.0  # [s]

    print('\t'.join(['t', 't_traj'] + [str(i) for i in range(otg.degrees_of_freedom)]))

    # Generate the trajectory within the control loop
    first_output, times, out_list = None, [], []
    res = Result.Working
    while res == Result.Working or res == Result.Paused:
        res = otg.update(inp, out)

        if out.time > 1.8 and phase == 'start':
            phase = 'brake'
        if res == Result.Paused and phase == 'brake':
            phase = 'accel'
        if phase == 'accel' and otg.speed >= 1.0:
            phase = 'end'

        if phase == 'brake':
            otg.speed = max(otg.speed - otg.delta_time / speed_change_duration, 0.0)
        if phase == 'accel':
            otg.speed = min(otg.speed + otg.delta_time / speed_change_duration, 1.0)

        # The out.time parameter denotes the time on the trajectory,
        # which is not the same as the time in the control loop as soon as the speed is not 1.0.
        time = (times[-1] if times else 0.0) + otg.delta_time

        print('\t'.join([f'{time:0.3f}', f'{out.time:0.3f}'] + [f'{p:0.3f}' for p in out.new_position]))
        times.append(time)
        out_list.append(copy(out))

        out.pass_to_input(inp)

        if not first_output:
            first_output = copy(out)

    print(f'Calculation duration: {first_output.calculation_duration:0.1f} [µs]')
    print(f'Trajectory duration: {first_output.trajectory.duration:0.4f} [s]')

    # Plot the trajectory
    # from pathlib import Path
    # from plotter import Plotter

    # project_path = Path(__file__).parent.parent.absolute()
    # Plotter.plot_trajectory(project_path / 'examples' / '16_trajectory.pdf', otg, inp, out_list, plot_jerk=False, times=times)
