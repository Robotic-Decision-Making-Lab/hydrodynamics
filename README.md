# Hydrodynamics

The hydrodynamics library provides a collection of C++ data classes for hydrodynamic parameters
(i.e., inertia, added mass, coriolis, damping, etc.). The goal of this project is to abstract away
hydrodynamic equations and representations so that you can instead focus on writing underwater
robotics algorithms.

## Installation

The hydrodynamics library has been implemented as a ROS 2 package with support for ROS 2 Humble,
Iron, and Rolling. To install and use this library in your own ROS 2 project, simply clone this
repository to your ROS 2 workspace:

```bash
cd path/to/ws_ros/src
git clone git@github.com:Robotic-Decision-Making-Lab/hydrodynamics.git
```

## Getting Started

All hydrodynamic parameters implemented in this project have been defined using Fossen's equations
for hydrodynamics. The hydrodynamic parameters implemented include:

- Mass: rigid body inertia and added mass
- Coriolis: centripetal, coriolis, and added coriolis
- Damping: linear and quadratic damping
- Restoring forces: buoyancy and gravitational forces
- Current effects: fluid velocity

Each of the aforementioned parameters provide their own distinct data class for independent use
or can be managed altogether within the `HydrodynamicParameters` class. For further
information regarding the interfaces available, please refer to the [library header](https://github.com/Robotic-Decision-Making-Lab/hydrodynamics/blob/main/include/hydrodynamics.hpp).

## License

The hydrodynamics library has been released under the MIT license.
