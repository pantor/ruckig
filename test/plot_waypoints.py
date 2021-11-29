from copy import copy
from pathlib import Path
from sys import path

from plotter import Plotter
import numpy as np

path.insert(0, str(Path(__file__).parent.absolute().parent / 'build'))

from ruckig import InputParameter, OutputParameter, Result, Ruckig, Trajectory, ControlInterface


def walk_through_trajectory(otg, inp):
    out_list = []
    out = OutputParameter(otg.degrees_of_freedom, 10)

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
        # [0.5, 1.6, 0],
        # [1, 1, 1],
        # [1.2, 2, 2],
        # [0.99, 0.2, -0.3],
        # [1.3, 0.5, 0],
        # [1.3, 0.5, 0],
        # [1.3, 0.5, 0],
        # [1.3, 0.5, 0],
        # [0.5, 1.6, 0],
        # [1, 1, 1],
        # [1.2, 2, 2],
        # [0.99, 0.2, -0.3],
        # [0.2, 0.1, 0.3],
        # [0.3, 0.2, 0.4],
        # [0.4, 0.3, 0.5],
        # [0.4, 0.4, 0.6],
        # [0.5, 0.5, 0.7],
        # [0.99, 0.9, 0.6],
    ]

    # inp.interrupt_calculation_duration = 500 # [µs]

    # np.random.seed(42)
    # inp.intermediate_positions = np.random.normal(0, 1, size=(10, 3))

    # inp.intermediate_positions = np.random.normal(0.19, 0.0, size=(5, 3))
    # inp.intermediate_positions = np.cumsum(inp.intermediate_positions, axis=0)

    # inp.intermediate_positions = np.array([
    #     np.linspace(0.01, 0.99, 100),
    #     np.linspace(0.01, 1.01, 100),
    #     np.linspace(0.02, 1.02, 100)
    # ]).T

    # inp.intermediate_positions = [
    #     [0.1],
    #     [0.9]
    # ]

    otg = Ruckig(inp.degrees_of_freedom, 0.01, 10)

    # out = OutputParameter(3)
    # np.random.seed(44)
    # durations = 0.0
    # for i in range(1000):
    #     inp.intermediate_positions = np.random.uniform(-2, 2, size=((i % 10), 3))
    #     traj = Trajectory(inp.degrees_of_freedom)
    #     otg.update(inp, out)
    #     durations += out.calculation_duration
    # print(f'Trajectory durations: {durations/1000:0.4f} [s]')

    out_list, time_offsets, time_offset = walk_through_trajectory(otg, inp)

    print(f'{0}\tCalculation duration: {out_list[0].calculation_duration:0.1f} [µs]')
    print(f'\tTrajectory duration: {out_list[0].trajectory.duration:0.4f} [s]')

    i = 1
    while out_list[i].new_calculation:
        print(f'{i}\tCalculation duration: {out_list[i].calculation_duration:0.1f} [µs]')
        print(f'\tTrajectory duration: {i*otg.delta_time + out_list[i].trajectory.duration:0.4f} [s]')
        i += 1

    Plotter.plot_trajectory('otg_trajectory.png', otg, inp, out_list, plot_jerk=False, time_offsets=time_offsets)
