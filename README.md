<div align="center">
  <h1 align="center">Ruckig</h1>
  <h3 align="center">
    Online Trajectory Generation. Real-time. Jerk-constrained.
  </h3>
</div>
<p align="center">
  <a href="https://github.com/pantor/ruckig/actions">
    <img src="https://github.com/pantor/ruckig/workflows/CI/badge.svg" alt="CI">
  </a>

  <a href="https://github.com/pantor/frankx/issues">
    <img src="https://img.shields.io/github/issues/pantor/frankx.svg" alt="Issues">
  </a>

  <a href="https://github.com/pantor/frankx/releases">
    <img src="https://img.shields.io/github/v/release/pantor/frankx.svg?include_prereleases&sort=semver" alt="Releases">
  </a>

  <a href="https://github.com/pantor/frankx/blob/master/LICENSE">
    <img src="https://img.shields.io/badge/license-MIT-green.svg" alt="LGPL">
  </a>
</p>


## Installation

```bash
mkdir -p build
cd build
cmake -DBUILD_TYPE=Release ..
make
make install
```

## Development

Frankx is written in C++17. It is currently tested against following versions

- Eigen v3.3.9
- Catch2 v2.13.3 (only for testing)
- Reflexxes v1.2.7 (only for testing)
- Pybind11 v2.6.0 (only for prototyping)
