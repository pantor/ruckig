from copy import copy
from pathlib import Path
from sys import path

from plotter import Plotter

path.insert(0, str(Path(__file__).parent.absolute().parent / 'build'))

from ruckig import InputParameter, OutputParameter, Result, Ruckig


def walk_through_trajectory(otg, inp, waypoints):
    out_list = []
    out = OutputParameter(inp.degrees_of_freedom)

    waypoints.append((inp.target_position, inp.target_velocity, inp.target_acceleration))

    time_offset = 0.0
    time_offsets = []
    for waypoint in waypoints:
        res = Result.Working
        while res == Result.Working:
            inp.target_position, inp.target_velocity, inp.target_acceleration = waypoint
            res = otg.update(inp, out)

            out.pass_to_input(inp)
            time_offsets.append(copy(time_offset))
            out_list.append(copy(out))

        time_offset += out_list[-1].trajectory.duration

    return out_list, time_offsets, time_offset


if __name__ == '__main__':
    inp = InputParameter(3)
    inp.current_position = [0, 0, 0]
    inp.current_velocity = [0, 0, 0]
    inp.current_acceleration = [0, 0, 0]
    inp.target_position = [1, 1, 1]
    inp.target_velocity = [0, 0, 0]
    inp.target_acceleration = [0, 0, 0]
    inp.max_velocity = [100, 100, 100]
    inp.max_acceleration = [100, 100, 100]
    inp.max_jerk = [1, 1, 1]

    intermediate_waypoints = [(
        [0.7, 0.2, 0.1],
        [0.56, 0.1, 0.1],
        [-0.3, 0, 0],
    )]

    otg = Ruckig(inp.degrees_of_freedom, 0.005)

    out_list, time_offsets, duration = walk_through_trajectory(otg, inp, intermediate_waypoints)


    print(f'Calculation duration: {out_list[0].calculation_duration:0.1f} [µs]')
    print(f'Trajectory duration: {duration:0.6f} [s]')

    Plotter.plot_trajectory('otg_multiple.png', otg, inp, out_list, plot_jerk=False, time_offsets=time_offsets)
