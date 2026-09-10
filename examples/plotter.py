from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


class Plotter:
    @staticmethod
    def _plot_data(  # ruff:ignore[too-many-arguments, too-many-branches, too-many-positional-arguments]
        path,
        inp,
        taxis,
        positions,
        velocities,
        accelerations,
        jerks=None,
        intermediate_durations=None,
        plot_acceleration=True,
        plot_jerk=False,
        title=None,
        show=False,
    ):
        plt.figure(figsize=(8.0, 2.0 + 3.0 * inp.degrees_of_freedom), dpi=120)
        plt.subplot(inp.degrees_of_freedom, 1, 1)
        if title:
            plt.title(title)

        for dof in range(inp.degrees_of_freedom):
            arrays = [positions[:, dof], velocities[:, dof], accelerations[:, dof]]
            if plot_jerk and jerks is not None:
                arrays.append(jerks[:, dof])
            global_max = np.max(arrays)
            global_min = np.min(arrays)

            plt.subplot(inp.degrees_of_freedom, 1, dof + 1)
            plt.ylabel(f'DoF {dof + 1}')
            plt.plot(taxis, positions[:, dof], label=f'Position {dof + 1}')
            plt.plot(taxis, velocities[:, dof], label=f'Velocity {dof + 1}')
            if plot_acceleration:
                plt.plot(taxis, accelerations[:, dof], label=f'Acceleration {dof + 1}')
            if plot_jerk and jerks is not None:
                plt.plot(taxis, jerks[:, dof], label=f'Jerk {dof + 1}')

            # Plot sections
            if intermediate_durations is not None:
                linewidth = 1.0 if len(intermediate_durations) < 20 else 0.25
                for t in intermediate_durations:
                    plt.axvline(x=t, color='black', linestyle='--', linewidth=linewidth)

            # Plot limit lines
            if inp.min_position[dof] > 1.4 * global_min:
                plt.axhline(y=inp.min_position[dof], color='tab:blue', linestyle='--', linewidth=1.1)
            if inp.max_position[dof] < 1.4 * global_max:
                plt.axhline(y=inp.max_position[dof], color='tab:blue', linestyle='--', linewidth=1.1)

            if inp.max_velocity[dof] < 1.4 * global_max:
                plt.axhline(y=inp.max_velocity[dof], color='tab:orange', linestyle='--', linewidth=1.1)
            min_velocity = inp.min_velocity[dof] if inp.min_velocity else -inp.max_velocity[dof]
            if min_velocity > 1.4 * global_min:
                plt.axhline(y=min_velocity, color='tab:orange', linestyle='--', linewidth=1.1)

            if plot_acceleration and inp.max_acceleration[dof] < 1.4 * global_max:
                plt.axhline(y=inp.max_acceleration[dof], color='tab:green', linestyle='--', linewidth=1.1)
            min_acceleration = inp.min_acceleration[dof] if inp.min_acceleration else -inp.max_acceleration[dof]
            if plot_acceleration and min_acceleration > 1.4 * global_min:
                plt.axhline(y=min_acceleration, color='tab:green', linestyle='--', linewidth=1.1)

            if plot_jerk and jerks is not None:
                if inp.max_jerk[dof] < 1.4 * global_max:
                    plt.axhline(y=inp.max_jerk[dof], color='tab:red', linestyle='--', linewidth=1.1)
                if -inp.max_jerk[dof] > 1.4 * global_min:
                    plt.axhline(y=-inp.max_jerk[dof], color='tab:red', linestyle='--', linewidth=1.1)

            plt.legend()
            plt.grid(True)

        plt.xlabel('t')
        plt.savefig(path)
        if show:
            plt.show()

    @staticmethod
    def plot(path: Path, trajectory, inp, show=False):
        taxis = np.linspace(0.0, trajectory.duration, num=500)
        positions, velocities, accelerations = [], [], []
        for t in taxis:
            position, velocity, acceleration = trajectory.at_time(t)
            positions.append(position)
            velocities.append(velocity)
            accelerations.append(acceleration)

        Plotter._plot_data(
            path, inp, taxis,
            np.array(positions), np.array(velocities), np.array(accelerations),
            intermediate_durations=trajectory.intermediate_durations,
            title=f'Trajectory with duration {trajectory.duration:.3f} s',
            show=show,
        )

    @staticmethod
    def plot_trajectory(  # ruff:ignore[too-many-arguments, too-many-positional-arguments]
        filename,
        ruckig,
        inp,
        out_list,
        show=False,
        plot_acceleration=True,
        plot_jerk=True,
        time_offsets=None,
        title=None,
        times=None,
    ):
        taxis = np.array(times if times else [x.time for x in out_list])
        if time_offsets:
            taxis += np.array(time_offsets)
        positions = np.array(list(map(lambda x: x.new_position, out_list)))
        velocities = np.array(list(map(lambda x: x.new_velocity, out_list)))
        accelerations = np.array(list(map(lambda x: x.new_acceleration, out_list)))

        # Numeric derivative to get jerk
        jerks = np.diff(accelerations, axis=0, prepend=accelerations[0, 0]) / ruckig.delta_time
        jerks[0, :] = 0.0
        jerks[-1, :] = 0.0

        intermediate_durations = None
        if hasattr(out_list[-1], 'trajectory'):
            intermediate_durations = out_list[-1].trajectory.intermediate_durations

        Plotter._plot_data(
            Path(__file__).parent.parent / 'build' / filename,
            inp, taxis, positions, velocities, accelerations,
            jerks=jerks,
            intermediate_durations=intermediate_durations,
            plot_acceleration=plot_acceleration,
            plot_jerk=plot_jerk,
            title=title,
            show=show,
        )
