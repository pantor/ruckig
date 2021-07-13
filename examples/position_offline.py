import time
from pathlib import Path
from sys import path

# Path to the build directory including a file similar to 'ruckig.cpython-37m-x86_64-linux-gnu'.
build_path = Path(__file__).parent.parent / 'build'
path.insert(0, str(build_path))

from ruckig import InputParameter, OutputParameter, Ruckig


if __name__ == '__main__':
    inp = InputParameter(3)
    inp.current_position = [0.2, 0, -1]
    inp.current_velocity = [0, 0.2, 0]
    inp.current_acceleration = [0, 1, 0]
    inp.target_position = [0, -1, -1]
    inp.target_velocity = [0.2, 0, 0]
    inp.target_acceleration = [0, 0.1, -0.1]
    inp.max_velocity = [2, 1, 1]
    inp.max_acceleration = [0.2, 2, 2]
    inp.max_jerk = [3, 4, 5]

    otg = Ruckig(3)

    out = OutputParameter(3)

    print('\t'.join(['t'] + [str(i) for i in range(otg.degrees_of_freedom)]))

    trajectory = otg.calculate(inp)
    print(f'Calculation duration: {out.calculation_duration:0.1f} [µs]')
    print(f'Trajectory duration: {out.trajectory.duration:0.4f} [s]')

    pos_updated = False
    t_start = time.time()
    t_now = time.time()-t_start
    while t_now < out.trajectory.duration:
        t_now = time.time()-t_start
        new_position, new_velocity, new_acceleration = out.trajectory.at_time(t_now)

        print('\t'.join([f'{t_now:0.3f}'] + [f'{p:0.3f}' for p in new_position]))

        if new_position[1] < -0.7 and not pos_updated:  # change target values on the fly and update trajectory
            pos_updated = True
            inp.current_position = new_position
            inp.current_velocity = new_velocity
            inp.current_acceleration = new_acceleration
            inp.target_position = [1, 0, 0]
            inp.max_velocity = [1, 2, 2]
            t_start = time.time()
            res = otg.update(inp, out)
            print(f'Calculation duration: {out.calculation_duration:0.1f} [µs]')
            print(f'Trajectory duration: {out.trajectory.duration:0.4f} [s]')

        time.sleep(0.02)  # limit output and calculation rate

    # verify that totally consumed time roughly matches trajectory.duration
    print(f'total time elapsed: {(time.time() - t_start)}')