// Copyright 2024, Evan Palmer
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#pragma once

#include <Eigen/Dense>
#include <expected>

namespace Eigen
{

// Extend the Eigen namespace to include commonly used matrix types
using Matrix6d = Eigen::Matrix<double, 6, 6>;
using Vector6d = Eigen::Matrix<double, 6, 1>;

}  // namespace Eigen

namespace hydrodynamics
{

struct Inertia
{
  Inertia() = default;

  /// Create a new wrapper for the inertial parameters (including added mass and rigid body inertia) using:
  /// - the mass of the vehicle,
  /// - the inertia tensor,
  /// - the added mass matrix,
  /// - the center of gravity of the vehicle.
  Inertia(double mass, const Eigen::Matrix3d & I, const Eigen::Matrix6d & added_mass, const Eigen::Vector3d & cog);

  auto operator*(const Eigen::Vector6d & accel) const -> Eigen::Vector6d { return mass_matrix * accel; }

  Eigen::Matrix6d rigid_body_matrix;
  Eigen::Matrix6d added_mass_matrix;
  Eigen::Matrix6d mass_matrix;
};

struct Coriolis
{
  Coriolis() = default;

  /// Create a new wrapper for the Coriolis and centripetal force parameters using:
  /// - the mass of the vehicle,
  /// - the inertia tensor,
  /// - the added mass matrix,
  /// - the center of gravity of the vehicle.
  Coriolis(double mass, const Eigen::Matrix3d & I, const Eigen::Matrix6d & added_mass, const Eigen::Vector3d & cog);

  /// Calculate the rigid body Coriolis matrix using the velocity of the vehicle.
  [[nodiscard]] auto rigid_body_coriolis_matrix(const Eigen::Vector6d & velocity) const -> Eigen::Matrix6d;

  /// Calculate the added mass Coriolis matrix using the velocity of the vehicle.
  [[nodiscard]] auto added_coriolis_matrix(const Eigen::Vector6d & velocity) const -> Eigen::Matrix6d;

  /// Calculate the Coriolis matrix using the velocity of the vehicle.
  [[nodiscard]] auto coriolis_matrix(const Eigen::Vector6d & velocity) const -> Eigen::Matrix6d;

  [[nodiscard]] auto operator()(const Eigen::Vector6d & velocity) const -> Eigen::Matrix6d
  {
    return coriolis_matrix(velocity);
  }

  Inertia inertia;
};

struct Damping
{
  Damping() = default;

  /// Create a new wrapper for the damping coefficients using:
  /// - the linear damping matrix,
  /// - the quadratic damping matrix.
  Damping(const Eigen::Matrix6d & linear, const Eigen::Matrix6d & quadratic);

  /// Create a new wrapper for the damping coefficients using:
  /// - the linear damping coefficients in the surge, sway, heave, roll, pitch, and yaw directions,
  /// - the quadratic damping coefficients in the surge, sway, heave, roll, pitch, and yaw directions.
  Damping(Eigen::Vector6d linear, Eigen::Vector6d quadratic);

  /// Calculate the damping matrix using the velocity of the vehicle.
  [[nodiscard]] auto damping_matrix(const Eigen::Vector6d & velocity) const -> Eigen::Matrix6d;

  [[nodiscard]] auto operator()(const Eigen::Vector6d & velocity) const -> Eigen::Matrix6d
  {
    return damping_matrix(velocity);
  }

  Eigen::Matrix6d linear_drag;
  Eigen::Matrix6d quadratic_drag;
};

struct RestoringForces
{
  RestoringForces() = default;

  /// Create a new wrapper for the restoring forces acting on the vehicle using:
  /// - the weight of the vehicle,
  /// - the buoyancy of the vehicle,
  /// - the center of buoyancy of the vehicle,
  /// - the center of gravity of the vehicle.
  RestoringForces(double weight, double buoyancy, Eigen::Vector3d cob, Eigen::Vector3d cog);

  /// Calculate the restoring forces acting on the vehicle using the current rotation of the vehicle.
  [[nodiscard]] auto restoring_forces_vector(const Eigen::Matrix3d & rotation) const -> Eigen::Vector6d;

  /// Calculate the restoring forces acting on the vehicle using the current rotation of the vehicle.
  [[nodiscard]] auto restoring_forces_vector(const Eigen::Quaterniond & rotation) const -> Eigen::Vector6d;

  [[nodiscard]] auto operator()(const Eigen::Matrix3d & rotation) const -> Eigen::Vector6d
  {
    return restoring_forces_vector(rotation);
  }

  double weight;
  double buoyancy;
  Eigen::Vector3d center_of_buoyancy;
  Eigen::Vector3d center_of_gravity;
};

struct Params
{
  Params() = default;

  /// Create a new wrapper for the hydrodynamic parameters using
  /// - the inertia of the vehicle,
  /// - the Coriolis and centripetal forces acting on the vehicle,
  /// - the damping coefficients acting on the vehicle,
  /// - the restoring forces acting on the vehicle.
  Params(const Inertia & M, const Coriolis & C, const Damping & D, const RestoringForces & g);

  /// Create a new wrapper for the hydrodynamic parameters using
  /// - the inertia of the vehicle,
  /// - the Coriolis and centripetal forces acting on the vehicle,
  /// - the damping coefficients acting on the vehicle,
  /// - the restoring forces acting on the vehicle.
  Params(Inertia && M, Coriolis && C, Damping && D, RestoringForces && g);

  /// Create a new wrapper for the hydrodynamic parameters using
  /// - the mass of the vehicle,
  /// - the inertia tensor,
  /// - the added mass matrix,
  /// - the linear damping matrix,
  /// - the quadratic damping matrix,
  /// - the center of gravity of the vehicle,
  /// - the center of buoyancy of the vehicle,
  /// - the weight of the vehicle,
  /// - the buoyancy of the vehicle.
  Params(
    double mass,
    const Eigen::Matrix3d & I,
    const Eigen::Matrix6d & added_mass,
    const Eigen::Matrix6d & linear_damping,
    const Eigen::Matrix6d & quadratic_damping,
    const Eigen::Vector3d & cog,
    const Eigen::Vector3d & cob,
    double weight,
    double buoyancy);

  Inertia M;
  Coriolis C;
  Damping D;
  RestoringForces g;
};

/// Perform forward dynamics to compute the acceleration of a rigid body using
/// - the hydrodynamic parameters of the vehicle,
/// - the vector of forces/torques acting on the vehicle (can also be a control input),
/// - the velocity of the vehicle,
/// - the orientation of the vehicle.
[[nodiscard]] auto forward_dynamics(
  const Params & params,
  const Eigen::Vector6d & vel,
  const Eigen::Vector6d & tau,
  const Eigen::Matrix3d & R) -> Eigen::Vector6d;

/// Perform inverse dynamics to compute the forces/torques acting on a rigid body using
/// - the hydrodynamic parameters of the vehicle,
/// - the acceleration of the vehicle,
/// - the velocity of the vehicle,
/// - the orientation of the vehicle.
[[nodiscard]] auto inverse_dynamics(
  const Params & params,
  const Eigen::Vector6d & acc,
  const Eigen::Vector6d & vel,
  const Eigen::Matrix3d & R) -> Eigen::Vector6d;

/// Parse the hydrodynamic model from a URDF string.
[[nodiscard]] auto parse_model_from_xml(const std::string & tree) -> std::expected<Params, std::string>;

/// Parse the hydrodynamic model from a URDF file.
[[nodiscard]] auto parse_model_from_urdf(const std::string & file) -> std::expected<Params, std::string>;

}  // namespace hydrodynamics
