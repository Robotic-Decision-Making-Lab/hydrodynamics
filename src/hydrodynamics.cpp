// Copyright 2023, Evan Palmer
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

#include "hydrodynamics/hydrodynamics.hpp"

#include <Eigen/Dense>

namespace hydrodynamics
{

namespace
{

auto make_skew_symmetric_matrix(const Eigen::Vector3d & v)
{
  Eigen::Matrix3d skew_symmetric_matrix;
  skew_symmetric_matrix << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
  return skew_symmetric_matrix;
}

}  // namespace

Inertia::Inertia(
  double mass,
  double Ixx,
  double Iyy,
  double Izz,
  double Xdu,
  double Ydv,
  double Zdw,
  double Kdp,
  double Mdq,
  double Ndr)
{
  // Construct the rigid body inertia matrix
  rigid_body_matrix = Eigen::Matrix6d::Zero();
  rigid_body_matrix.topLeftCorner(3, 3) = mass * Eigen::Matrix3d::Identity();
  rigid_body_matrix.bottomRightCorner(3, 3) = Eigen::Vector3d(Ixx, Iyy, Izz).asDiagonal().toDenseMatrix();

  // Construct the added mass matrix
  added_mass_matrix = -Eigen::Vector6d(Xdu, Ydv, Zdw, Kdp, Mdq, Ndr).asDiagonal().toDenseMatrix();

  // The complete mass matrix is the sum of the rigid body and added mass matrices
  mass_matrix = rigid_body_matrix + added_mass_matrix;
}

Inertia::Inertia(double mass, const Eigen::Vector3d & moments, const Eigen::Vector6d & added_mass)
{
  // Construct the rigid body inertia matrix
  Eigen::Matrix6d rigid_body_matrix = Eigen::Matrix6d::Zero();
  rigid_body_matrix.topLeftCorner(3, 3) = mass * Eigen::Matrix3d::Identity();
  rigid_body_matrix.bottomRightCorner(3, 3) = moments.asDiagonal().toDenseMatrix();

  // Construct the added mass matrix
  added_mass_matrix = -added_mass.asDiagonal().toDenseMatrix();

  // The complete mass matrix is the sum of the rigid body and added mass matrices
  mass_matrix = rigid_body_matrix + added_mass_matrix;
}

Coriolis::Coriolis(
  double mass,
  double Ixx,
  double Iyy,
  double Izz,
  double Xdu,
  double Ydv,
  double Zdw,
  double Kdp,
  double Mdq,
  double Ndr)
: mass(mass),
  moments(Eigen::Vector3d(Ixx, Iyy, Izz).asDiagonal().toDenseMatrix()),
  added_mass_coeff(Eigen::Vector6d(Xdu, Ydv, Zdw, Kdp, Mdq, Ndr))
{
}

Coriolis::Coriolis(double mass, const Eigen::Vector3d & moments, Eigen::Vector6d added_mass)
: mass(mass),
  moments(moments.asDiagonal().toDenseMatrix()),
  added_mass_coeff(std::move(added_mass))
{
}

auto Coriolis::calculate_rigid_body_coriolis_matrix(const Eigen::Vector3d & angular_velocity) const -> Eigen::Matrix6d
{
  Eigen::Matrix6d coriolis = Eigen::Matrix6d::Zero();

  const Eigen::Vector3d moments_v2 = moments * angular_velocity;

  coriolis.topLeftCorner(3, 3) = mass * make_skew_symmetric_matrix(angular_velocity);
  coriolis.bottomRightCorner(3, 3) = -make_skew_symmetric_matrix(moments_v2);

  return coriolis;
}

auto Coriolis::calculate_added_coriolis_matrix(const Eigen::Vector6d & velocity) const -> Eigen::Matrix6d
{
  Eigen::Matrix6d coriolis = Eigen::Matrix6d::Zero();

  const Eigen::Matrix3d linear_vel = make_skew_symmetric_matrix(Eigen::Vector3d(
    added_mass_coeff(0) * velocity(0), added_mass_coeff(1) * velocity(1), added_mass_coeff(2) * velocity(2)));

  const Eigen::Matrix3d angular_vel = make_skew_symmetric_matrix(Eigen::Vector3d(
    added_mass_coeff(3) * velocity(3), added_mass_coeff(4) * velocity(4), added_mass_coeff(5) * velocity(5)));

  coriolis.topRightCorner(3, 3) = linear_vel;
  coriolis.bottomLeftCorner(3, 3) = linear_vel;
  coriolis.bottomRightCorner(3, 3) = angular_vel;

  return coriolis;
}

auto Coriolis::calculate_coriolis_matrix(const Eigen::Vector6d & velocity) const -> Eigen::Matrix6d
{
  return calculate_rigid_body_coriolis_matrix(velocity.bottomRows(3)) + calculate_added_coriolis_matrix(velocity);
}

Damping::Damping(
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
  double Nrr)
: linear_coeff(Eigen::Vector6d(Xu, Yv, Zw, Kp, Mq, Nr)),
  quadratic_coeff(Eigen::Vector6d(Xuu, Yvv, Zww, Kpp, Mqq, Nrr))
{
}

Damping::Damping(Eigen::Vector6d linear, Eigen::Vector6d quadratic)
: linear_coeff(std::move(linear)),
  quadratic_coeff(std::move(quadratic))
{
}

auto Damping::calculate_damping_matrix(const Eigen::Vector6d & velocity) const -> Eigen::Matrix6d
{
  const Eigen::Matrix6d quadratic =
    -(quadratic_coeff.asDiagonal().toDenseMatrix() * velocity.cwiseAbs()).asDiagonal().toDenseMatrix();

  return -linear_coeff.asDiagonal().toDenseMatrix() + quadratic;
}

RestoringForces::RestoringForces(
  double weight,
  double buoyancy,
  Eigen::Vector3d center_of_buoyancy,
  Eigen::Vector3d center_of_gravity)
: weight(weight),
  buoyancy(buoyancy),
  center_of_buoyancy(std::move(center_of_buoyancy)),
  center_of_gravity(std::move(center_of_gravity))
{
}

auto RestoringForces::calculate_restoring_forces_vector(const Eigen::Matrix3d & rotation) const -> Eigen::Vector6d
{
  const Eigen::Vector3d fg(0, 0, weight);
  const Eigen::Vector3d fb(0, 0, -buoyancy);

  Eigen::Vector6d g_rb;

  const Eigen::Matrix3d rot_t = rotation.transpose();

  g_rb.topRows(3) = rot_t * (fg + fb);
  g_rb.bottomRows(3) = center_of_gravity.cross(rot_t * fg) + center_of_buoyancy.cross(rot_t * fb);

  g_rb *= -1;

  return g_rb;
}

}  // namespace hydrodynamics
