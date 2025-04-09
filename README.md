# Hydrodynamics

The hydrodynamics library provides a high-level API for storing/parsing
hydrodynamic parameters and performing forward/inverse dynamics.

## Installation

hydrodynamics supports ROS 2 Jazzy, and Rolling. To install and use
this library in your own project, clone this repository to your workspace:

```bash
cd path/to/ws_ros/src
git clone git@github.com:Robotic-Decision-Making-Lab/hydrodynamics.git
```

## Getting Started

All parameters implemented in this project have been defined using Fossen's
equations for hydrodynamics, i.e.,

```math
M\dot{v} + C(v)v + D(v)v + g(\eta) = \tau,
```

where $M$ is the mass matrix, including rigid body inertia and added mass;
$C(v)$ is the matrix of Coriolis and centripetal effects; $D(v)$ is the drag
matrix, including linear and quadratic drag; and $g(\eta)$ is the vector of
restoring forces.

For more information regarding how to use the library, please see the project
[examples](https://github.com/Robotic-Decision-Making-Lab/hydrodynamics/tree/main/examples).

## License

The hydrodynamics library has been released under the MIT license.
