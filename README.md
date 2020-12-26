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

```bash
mkdir -p build
cd build
cmake -DBUILD_TYPE=Release ..
make
make install
```

A python module can be built using the `BUILD_PYTHON_MODULE` CMake flag.


## Tutorial

Figure.


### Real-time trajectory generation

```c++
Ruckig<6> ruckig {0.001}; // DoF, control cycle in [s]

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

### Input Type


### Result Type


### Output Type


### Exceptions


## Development

Ruckig is written in C++17. It is currently tested against following versions

- Eigen v3.3.9
- Catch2 v2.13.3 (only for testing)
- Reflexxes v1.2.7 (only for testing)
- Pybind11 v2.6.0 (only for prototyping)

The current test suite validates over 170.000 trajectories based on random inputs.
