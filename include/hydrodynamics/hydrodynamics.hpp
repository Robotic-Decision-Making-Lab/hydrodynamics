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

#include "hydrodynamics/eigen.hpp"

namespace hydrodynamics
{

struct Inertia
{
  Inertia() = default;

  /// Create a new wrapper for inertial parameters (including added mass and rigid body inertia) using:
  /// - the mass of the vehicle,
  /// - the moments of inertia about the x, y, and z axes,
  /// - the added mass coefficients in the surge, sway, heave, roll, pitch, and yaw directions.
  Inertia(
    double mass,
    double Ixx,
    double Iyy,
    double Izz,
    double Xdu,
    double Ydv,
    double Zdw,
    double Kdp,
    double Mdq,
    double Ndr);

  /// Create a new wrapper for inertial parameters (including added mass and rigid body inertia) using:
  /// - the mass of the vehicle,
  /// - the moments of inertia about the x, y, and z axes,
  /// - the added mass coefficients in the surge, sway, heave, roll, pitch, and yaw directions.
  Inertia(double mass, const Eigen::Vector3d & moments, const Eigen::Vector6d & added_mass);

  Eigen::Matrix6d rigid_body_matrix;
  Eigen::Matrix6d added_mass_matrix;
  Eigen::Matrix6d mass_matrix;
};

struct Coriolis
{
  Coriolis() = default;

  /// Create a new wrapper for the Coriolis and centripetal force parameters using:
  /// - the mass of the vehicle,
  /// - the moments of inertia about the x, y, and z axes,
  /// - the added mass coefficients in the surge, sway, heave, roll, pitch, and yaw directions.
  Coriolis(
    double mass,
    double Ixx,
    double Iyy,
    double Izz,
    double Xdu,
    double Ydv,
    double Zdw,
    double Kdp,
    double Mdq,
    double Ndr);

  /// Create a new wrapper for the Coriolis and centripetal force parameters using:
  /// - the mass of the vehicle,
  /// - the moments of inertia about the x, y, and z axes,
  /// - the added mass coefficients in the surge, sway, heave, roll, pitch, and yaw directions.
  Coriolis(double mass, const Eigen::Vector3d & moments, Eigen::Vector6d added_mass);

  /// Calculate the rigid body Coriolis matrix using the angular velocity of the vehicle.
  [[nodiscard]] auto calculate_rigid_body_coriolis_matrix(const Eigen::Vector3d & angular_velocity) const
    -> Eigen::Matrix6d;

  /// Calculate the added mass Coriolis matrix using the velocity of the vehicle.
  [[nodiscard]] auto calculate_added_coriolis_matrix(const Eigen::Vector6d & velocity) const -> Eigen::Matrix6d;

  /// Calculate the Coriolis matrix using the velocity of the vehicle.
  [[nodiscard]] auto calculate_coriolis_matrix(const Eigen::Vector6d & velocity) const -> Eigen::Matrix6d;

  double mass;
  Eigen::Matrix3d moments;
  Eigen::Vector6d added_mass_coeff;
};

struct Damping
{
  Damping() = default;

  /// Create a new wrapper for the damping coefficients using:
  /// - the linear damping coefficients in the surge, sway, heave, roll, pitch, and yaw directions,
  /// - the quadratic damping coefficients in the surge, sway, heave, roll, pitch, and yaw directions.
  Damping(
    double Xu,
    double Yv,
    double Zw,
    double Kp,
    double Mq,
    double Nr,
    double Xuu,
    double Yvv,
    double Zww,
    double Kpp,
    double Mqq,
    double Nrr);

  /// Create a new wrapper for the damping coefficients using:
  /// - the linear damping coefficients in the surge, sway, heave, roll, pitch, and yaw directions,
  /// - the quadratic damping coefficients in the surge, sway, heave, roll, pitch, and yaw directions.
  Damping(Eigen::Vector6d linear, Eigen::Vector6d quadratic);

  /// Calculate the damping matrix using the velocity of the vehicle.
  [[nodiscard]] auto calculate_damping_matrix(const Eigen::Vector6d & velocity) const -> Eigen::Matrix6d;

  Eigen::Vector6d linear_coeff;
  Eigen::Vector6d quadratic_coeff;
};

struct RestoringForces
{
  RestoringForces() = default;

  /// Create a new wrapper for the restoring forces acting on the vehicle using:
  /// - the weight of the vehicle,
  /// - the buoyancy of the vehicle,
  /// - the center of buoyancy of the vehicle,
  /// - the center of gravity of the vehicle.
  RestoringForces(
    double weight,
    double buoyancy,
    Eigen::Vector3d center_of_buoyancy,
    Eigen::Vector3d center_of_gravity);

  /// Calculate the restoring forces acting on the vehicle using the current rotation of the vehicle.
  [[nodiscard]] auto calculate_restoring_forces_vector(const Eigen::Matrix3d & rotation) const -> Eigen::Vector6d;

  double weight;
  double buoyancy;
  Eigen::Vector3d center_of_buoyancy;
  Eigen::Vector3d center_of_gravity;
};

}  // namespace hydrodynamics
