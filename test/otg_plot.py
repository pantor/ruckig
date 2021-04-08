from copy import copy
from pathlib import Path
from sys import path

import matplotlib.pyplot as plt
import numpy as np

path.insert(0, str(Path(__file__).parent.parent / 'build'))

from ruckig import Quintic, InputParameter, OutputParameter, Result, Ruckig, Smoothie, Synchronization, Interface, DurationDiscretization
from ruckig import Reflexxes


def walk_through_trajectory(otg, inp, print_table=True):
    t_list, out_list = [], []
    out = OutputParameter(3)

    res = Result.Working
    old_acc = 0
    print_dof = 0
    while res == Result.Working:
        res = otg.update(inp, out)

        inp.current_position = out.new_position
        inp.current_velocity = out.new_velocity
        inp.current_acceleration = out.new_acceleration

        if print_table:
            jerk = (old_acc - out.new_acceleration[print_dof]) / otg.delta_time
            old_acc = out.new_acceleration[print_dof]
            # print(str(out.time) + '\t' + str(inp.current_position[print_dof]) + '\t' + str(inp.current_velocity[print_dof]) + '\t' + str(inp.current_acceleration[print_dof]) + '\t' + str(jerk))
            # print(str(inp.current_position[0]) + '\t' + str(inp.current_position[1]))

        t_list.append(out.time)
        out_list.append(copy(out))

    return t_list, out_list


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
    # inp.interface = Interface.Velocity
    # inp.synchronization = Synchronization.Phase
    # inp.duration_discretization = DurationDiscretization.Discrete

    inp.current_position = [-4.490717417930574, 3.467236624628543, -0.7545929089757601]
    inp.current_velocity = [0.1839756723363622, -0.4356283320280516, 0.7490399525818022]
    inp.current_acceleration = [-1.057769973808928, 0, -2.368645439140517]
    inp.target_position = [-4.928244836531066, -4.821780824003112, -8.20567952461017]
    inp.target_velocity = [0.1097319156272965, -0.9272874846270881, 0]
    inp.target_acceleration = [0.03089046366221739, -0.9744054582899561, 0]
    inp.max_velocity = [6.144314006624488, 2.93258338415229, 0.1820021269527196]
    inp.max_acceleration = [5.199401036221791, 1.848176490768948, 11.11168017805234]
    inp.max_jerk = [9.940940357283978, 10.46997753899755, 0.08166297169205029]

    # otg = Quintic(3, 0.005)
    # otg = Smoothie(3, 0.005)
    # otg = Reflexxes(3, 0.005)
    otg = Ruckig(3, 0.5)

    t_list, out_list = walk_through_trajectory(otg, inp)

    # print(out_list[0].trajectory.get_position_extrema())

    print(f'Calculation duration: {out_list[0].calculation_duration:0.1f} [µs]')
    print(f'Trajectory duration: {out_list[0].trajectory.duration:0.4f} [s]')

    plot_trajectory(t_list, out_list)
