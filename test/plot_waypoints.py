from copy import copy
from pathlib import Path
from sys import path

from plotter import Plotter
import numpy as np

path.insert(0, str(Path(__file__).parent.absolute().parent / 'build'))

from ruckig import InputParameter, OutputParameter, Result, Ruckig, Trajectory, ControlInterface


def walk_through_trajectory(otg, inp):
    out_list = []
    out = OutputParameter(otg.degrees_of_freedom, otg.max_number_of_waypoints)

    time_offset = 0.0
    time_offsets = []

    res = Result.Working
    while res == Result.Working:
        res = otg.update(inp, out)
        out_list.append(copy(out))
        out.pass_to_input(inp)

        time_offset += out.time if out.new_calculation else 0.0
        time_offsets.append(copy(time_offset))

    return out_list, time_offsets, time_offset


if __name__ == '__main__':
    inp = InputParameter(3)

    inp.current_position = [0, 0, 0]
    inp.current_velocity = [0.1, 0, 0]
    inp.current_acceleration = [0, 0, 0]
    inp.target_position = [1, 1, 1]
    inp.target_velocity = [0, 0, 0]
    inp.target_acceleration = [0, 0, 0]
    inp.max_velocity = [2, 2, 2]
    inp.max_acceleration = [2, 2, 2]
    inp.max_jerk = [2, 2, 2]

    # inp.control_interface = ControlInterface.Velocity
    # inp.minimum_duration = 7.0
    # inp.enabled = [False, True, False]

    inp.intermediate_positions = [
        [-1, -1, 0],
        [1.3, 0.5, 0],
    ]

    otg = Ruckig(inp.degrees_of_freedom, 0.01, 10)

    # otg.calculator.waypoints_calculator.number_global_steps = 1
    # otg.calculator.waypoints_calculator.number_local_steps = 16
    # otg.calculator.waypoints_calculator.number_smoothing_steps = 0
    # otg.calculator.waypoints_calculator.number_acceleration_smoothing_steps = 0

    out_list, time_offsets, time_offset = walk_through_trajectory(otg, inp)

    print(f'{0}\tCalculation duration: {out_list[0].calculation_duration:0.1f} [µs]')
    print(f'\tTrajectory duration: {out_list[0].trajectory.duration:0.4f} [s]')

    i = 1
    while out_list[i].new_calculation:
        print(f'{i}\tCalculation duration: {out_list[i].calculation_duration:0.1f} [µs]')
        print(f'\tTrajectory duration: {i*otg.delta_time + out_list[i].trajectory.duration:0.4f} [s]')
        i += 1

    Plotter.plot_trajectory('otg_trajectory.png', otg, inp, out_list, plot_jerk=False, time_offsets=time_offsets)
