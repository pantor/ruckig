<div align="center">
  <h1 align="center">Ruckig</h1>
  <h3 align="center">
    Online Trajectory Generation. Real-time. Time-optimal. Jerk-constrained.
  </h3>
</div>
<p align="center">
  <a href="https://github.com/pantor/ruckig/actions">
    <img src="https://github.com/pantor/ruckig/workflows/CI/badge.svg" alt="CI">
  </a>

  <a href="https://github.com/pantor/ruckig/issues">
    <img src="https://img.shields.io/github/issues/pantor/ruckig.svg" alt="Issues">
  </a>

  <a href="https://github.com/pantor/ruckig/releases">
    <img src="https://img.shields.io/github/v/release/pantor/ruckig.svg?include_prereleases&sort=semver" alt="Releases">
  </a>

  <a href="https://github.com/pantor/ruckig/blob/master/LICENSE">
    <img src="https://img.shields.io/badge/license-MIT-green.svg" alt="LGPL">
  </a>
</p>

Ruckig calculates a time-optimal trajectory given a *target* waypoint with position, velocity, and acceleration, starting from *any* initial state limited by velocity, acceleration, and jerk constraints. Robotics. Machine control. Ruckig is a more powerful and open-source alternative to the [Reflexxes Type IV](http://reflexxes.ws/) library. In fact, Ruckig is a Type V trajectory generator. In general, Ruckig allows for instant reactions to unforeseen events.


## Installation

To build Ruckig using CMake, just 

```bash
mkdir -p build
cd build
cmake -DBUILD_TYPE=Release ..
make
```

To install Ruckig in a system-wide directory, use `(sudo) make install`. We recommend to include

A python module can be built using the `BUILD_PYTHON_MODULE` CMake flag.


## Tutorial

Figure. Currently only the (more-complex) *position* interface is implemented.

![Trajectory Profile](/doc/example_profile.png?raw=true)

### Real-time trajectory generation

```c++
Ruckig<6> ruckig {0.001}; // Number DoFs; control cycle in [s]

InputParameter<6> input;
input.current_position = {};
input.current_velocity = {};
input.current_acceleration = {};
input.target_position = {};
input.target_velocity = {};
input.target_acceleration = {};
input.max_velocity = {};
input.max_acceleration = {};
input.max_jerk = {};

OutputParameter<6> output;

while (otg.update(input, output) == Result::Working) {
  // output.new_position

  input.current_position = output.new_position;
  input.current_velocity = output.new_velocity;
  input.current_acceleration = output.new_acceleration;
}

```

`at_time(double time)`


### Input Type

The input type 

```c++
std::array<double, DOFs> current_position;
std::array<double, DOFs> current_velocity; // Initialized to zero
std::array<double, DOFs> current_acceleration; // Initialized to zero

std::array<double, DOFs> target_position;
std::array<double, DOFs> target_velocity; // Initialized to zero
std::array<double, DOFs> target_acceleration; // Initialized to zero

std::array<double, DOFs> max_velocity;
std::array<double, DOFs> max_acceleration;
std::array<double, DOFs> max_jerk;

std::array<bool, DOFs> enabled; // Initialized to true
std::optional<double> minimum_duration;
```

To check the input in front of a calculation step, the `ruckig.validate_input(input)` method returns `false` if an input is not valid.


### Result Type

The `update` function of the Ruckig class returns a Result type that indicates the current state of the algorithm. Currently, this can either be **working**, **finished** if the trajectory has finished, or **error** if something went wrong during calculation. In this case, an exception (see below for more details) is thrown.

State    | Error Code
-------- | ----------
Working  | 0
Finished | 1
Error    | -1
ErrorInvalidInput           | -100
ErrorExecutionTime          | -101
ErrorSynchronization        | -102
ErrorNoPhaseSynchronization | -103
ErrorUserTimeOutOfRange     | -105
-------- | ------------


### Output Type

```c++
std::array<double, DOFs> new_position;
std::array<double, DOFs> new_velocity;
std::array<double, DOFs> new_acceleration;

double duration; // Duration of the trajectory [s]
bool new_calculation; // Whether a new calactuion was performed in the last cycle
double calculation_duration; // Duration of the calculation in the last cycle [µs]
```


## Tests

The current test suite validates over 190.000 (random) trajectories. The numerical exactness is tested for the position, velocity, acceleration, and time target to be within `1e-8`, for the velocity and acceleration limit to be withing `1e-9`, and for the jerk limit to be within a numerical error of `1e-12`.


## Development

Ruckig is written in C++17. It is currently tested against following versions

- Eigen v3.3.9
- Catch2 v2.13.3 (only for testing)
- Reflexxes v1.2.7 (only for testing)
- Pybind11 v2.6.0 (only for prototyping)
