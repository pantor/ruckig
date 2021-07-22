from copy import copy
from pathlib import Path
from sys import path

import matplotlib.pyplot as plt
import numpy as np

path.insert(0, str(Path(__file__).parent.parent / 'build'))

from ruckig import OutputParameter, Result


class Plotter:
    @staticmethod
    def plot_trajectory(filename, otg, inp, t_list, out_list, show=False):
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
        plt.savefig(Path(__file__).parent.parent / 'build' / filename)

        if show:
            plt.show()
