from copy import copy
from pathlib import Path
from sys import path

import matplotlib.pyplot as plt
import numpy as np

path.insert(0, str(Path(__file__).parent.parent / 'build'))

from ruckig import InputParameter, OutputParameter, Result, Ruckig, Synchronization, Interface, DurationDiscretization
from ruckig import Reflexxes


def walk_through_trajectory(otg, inp, intermediate_targets):
    t_list, out_list = [], []
    out = OutputParameter(3)

    old_target = inp.target_position, inp.target_velocity, inp.target_acceleration
    intermediate_targets.append(old_target)
    
    time_offset = 0.0
    for inp.target_position, inp.target_velocity, inp.target_acceleration in intermediate_targets:
        res = Result.Working
        while res == Result.Working:
            res = otg.update(inp, out)

            inp.current_position = out.new_position
            inp.current_velocity = out.new_velocity
            inp.current_acceleration = out.new_acceleration

            t_list.append(time_offset + out.time)
            out_list.append(copy(out))

        time_offset += out.trajectory.duration

    return t_list, out_list, time_offset


def plot_trajectory(t_list, out_list):
    qaxis = np.array(list(map(lambda x: x.new_position, out_list)))
    dqaxis = np.array(list(map(lambda x: x.new_velocity, out_list)))
    ddqaxis = np.array(list(map(lambda x: x.new_acceleration, out_list)))
    dddqaxis = np.diff(ddqaxis, axis=0, prepend=ddqaxis[0, 0]) / otg.delta_time
    dddqaxis[0, :] = 0.0
    dddqaxis[-1, :] = 0.0

    plt.figure(figsize=(8.0, 2.0 + 3.0 * inp.degrees_of_freedom), dpi=120)

    for dof in range(inp.degrees_of_freedom):
        global_max = np.max([qaxis[:, dof], dqaxis[:, dof], ddqaxis[:, dof], dddqaxis[:, dof]])
        global_min = np.min([qaxis[:, dof], dqaxis[:, dof], ddqaxis[:, dof], dddqaxis[:, dof]])

        plt.subplot(inp.degrees_of_freedom, 1, dof + 1)
        plt.plot(t_list, qaxis[:, dof], label=f'Position {dof+1}')
        plt.plot(t_list, dqaxis[:, dof], label=f'Velocity {dof+1}')
        plt.plot(t_list, ddqaxis[:, dof], label=f'Acceleration {dof+1}')
        plt.plot(t_list, dddqaxis[:, dof], label=f'Jerk {dof+1}')

        # Plot limit lines
        if inp.max_velocity[dof] < 1.4 * global_max:
            plt.axhline(y=inp.max_velocity[dof], color='orange', linestyle='--', linewidth=1.1)

        min_velocity = inp.min_velocity[dof] if inp.min_velocity else -inp.max_velocity[dof]
        if min_velocity > 1.4 * global_min:
            plt.axhline(y=min_velocity, color='orange', linestyle='--', linewidth=1.1)

        if inp.max_acceleration[dof] < 1.4 * global_max:
            plt.axhline(y=inp.max_acceleration[dof], color='g', linestyle='--', linewidth=1.1)

        min_acceleration = inp.min_acceleration[dof] if inp.min_acceleration else -inp.max_acceleration[dof]
        if min_acceleration > 1.4 * global_min:
            plt.axhline(y=min_acceleration, color='g', linestyle='--', linewidth=1.1)

        if inp.max_jerk[dof] < 1.4 * global_max:
            plt.axhline(y=inp.max_jerk[dof], color='red', linestyle='--', linewidth=1.1)

        if -inp.max_jerk[dof] > 1.4 * global_min:
            plt.axhline(y=-inp.max_jerk[dof], color='red', linestyle='--', linewidth=1.1)

        plt.legend()
        plt.grid(True)

    plt.xlabel('t')
    plt.savefig(Path(__file__).parent.parent / 'build' / 'otg_trajectory.png')
    # plt.show()


if __name__ == '__main__':
    inp = InputParameter(3)
    inp.current_position = [0, 0, 0]
    inp.current_velocity = [0, 0, 0]
    inp.current_acceleration = [0, 0, 0]
    inp.target_position = [1, 0, 0]
    inp.target_velocity = [0, 0, 0]
    inp.target_acceleration = [0, 0, 0]
    inp.max_velocity = [100, 100, 100]
    inp.max_acceleration = [100, 100, 100]
    inp.max_jerk = [1, 1, 1]

    intermediate_targets = [(
        [0.7, 0.0, 0.0],
        [0.56, 0.0, 0.0],
        [-0.3, 0, 0],
    )]

    otg = Ruckig(3, 0.005)

    t_list, out_list, duration = walk_through_trajectory(otg, inp, intermediate_targets)


    print(f'Calculation duration: {out_list[0].calculation_duration:0.1f} [µs]')
    print(f'Trajectory duration: {duration:0.4f} [s]')

    plot_trajectory(t_list, out_list)
