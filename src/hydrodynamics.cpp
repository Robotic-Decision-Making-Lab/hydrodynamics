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

#include "hydrodynamics.hpp"

#include <Eigen/Dense>

namespace hydrodynamics
{

namespace
{

Eigen::Matrix3d makeSkewSymmetricMatrix(const Eigen::Vector3d & v)
{
  Eigen::Matrix3d skew_symmetric_matrix;
  skew_symmetric_matrix << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
  return skew_symmetric_matrix;
}

}  // namespace

Inertia::Inertia(double mass, double Ixx, double Iyy, double Izz, double Xdu, double Ydv, double Zdw, double Kdp,
                 double Mdq, double Ndr)
{
  rigid_body_mat_ = Eigen::Matrix6d::Zero();
  rigid_body_mat_.topLeftCorner(3, 3) = mass * Eigen::Matrix3d::Identity();
  rigid_body_mat_.bottomRightCorner(3, 3) = Eigen::Vector3d(Ixx, Iyy, Izz).asDiagonal().toDenseMatrix();

  added_mass_mat_ = -Eigen::Vector6d(Xdu, Ydv, Zdw, Kdp, Mdq, Ndr).asDiagonal().toDenseMatrix();

  inertia_mat_ = rigid_body_mat_ + added_mass_mat_;
}

Inertia::Inertia(double mass, double Ixx, double Ixy, double Ixz, double Iyy, double Iyz, double Izz, double Xdu,
                 double Ydv, double Zdw, double Kdp, double Mdq, double Ndr)
{
  rigid_body_mat_ = Eigen::Matrix6d::Zero();
  rigid_body_mat_.topLeftCorner(3, 3) = mass * Eigen::Matrix3d::Identity();
  rigid_body_mat_.bottomRightCorner(3, 3) << Ixx, Ixy, Ixz, Ixy, Iyy, Iyz, Ixz, Iyz, Izz;

  added_mass_mat_ = -Eigen::Vector6d(Xdu, Ydv, Zdw, Kdp, Mdq, Ndr).asDiagonal().toDenseMatrix();

  inertia_mat_ = rigid_body_mat_ + added_mass_mat_;
}

Inertia::Inertia(double mass, Eigen::Vector3d moments, Eigen::Vector6d added_mass)
{
  rigid_body_mat_ = Eigen::Matrix6d::Zero();
  rigid_body_mat_.topLeftCorner(3, 3) = mass * Eigen::Matrix3d::Identity();
  rigid_body_mat_.bottomRightCorner(3, 3) = moments.asDiagonal().toDenseMatrix();

  added_mass_mat_ = -added_mass.asDiagonal().toDenseMatrix();

  inertia_mat_ = rigid_body_mat_ + added_mass_mat_;
}

Inertia::Inertia(double mass, Eigen::Matrix3d moments, Eigen::Vector6d added_mass)
{
  rigid_body_mat_ = Eigen::Matrix6d::Zero();
  rigid_body_mat_.topLeftCorner(3, 3) = mass * Eigen::Matrix3d::Identity();
  rigid_body_mat_.bottomRightCorner(3, 3) = moments;

  added_mass_mat_ = -added_mass.asDiagonal().toDenseMatrix();

  inertia_mat_ = rigid_body_mat_ + added_mass_mat_;
}

Coriolis::Coriolis(double mass, double Ixx, double Iyy, double Izz, double Xdu, double Ydv, double Zdw, double Kdp,
                   double Mdq, double Ndr)
: mass_(mass),
  moments_mat_(Eigen::Vector3d(Ixx, Iyy, Izz).asDiagonal().toDenseMatrix()),
  added_mass_coeff_(Eigen::Vector6d(Xdu, Ydv, Zdw, Kdp, Mdq, Ndr))
{
}

Coriolis::Coriolis(double mass, double Ixx, double Ixy, double Ixz, double Iyy, double Iyz, double Izz, double Xdu,
                   double Ydv, double Zdw, double Kdp, double Mdq, double Ndr)
: mass_(mass), added_mass_coeff_(Eigen::Vector6d(Xdu, Ydv, Zdw, Kdp, Mdq, Ndr))
{
  moments_mat_ = Eigen::Matrix3d::Zero();
  moments_mat_ << Ixx, Ixy, Ixz, Ixy, Iyy, Iyz, Ixz, Iyz, Izz;
}

Coriolis::Coriolis(double mass, Eigen::Vector3d moments, Eigen::Vector6d added_mass)
: mass_(mass), moments_mat_(moments.asDiagonal().toDenseMatrix()), added_mass_coeff_(std::move(added_mass))
{
}

Coriolis::Coriolis(double mass, Eigen::Matrix3d moments, Eigen::Vector6d added_mass)
: mass_(mass), moments_mat_(std::move(moments)), added_mass_coeff_(std::move(added_mass))
{
}

Eigen::Matrix6d Coriolis::calculateRigidBodyCoriolisMatrix(const Eigen::Vector3d & angular_velocity) const
{
  Eigen::Matrix6d coriolis = Eigen::Matrix6d::Zero();

  const Eigen::Vector3d moments_v2 = moments_mat_ * angular_velocity;

  coriolis.topLeftCorner(3, 3) = mass_ * makeSkewSymmetricMatrix(angular_velocity);
  coriolis.bottomRightCorner(3, 3) = -makeSkewSymmetricMatrix(moments_v2);

  return coriolis;
}

Eigen::Matrix6d Coriolis::calculateAddedCoriolisMatrix(const Eigen::Vector6d & velocity) const
{
  Eigen::Matrix6d coriolis = Eigen::Matrix6d::Zero();

  Eigen::Matrix3d linear_vel = makeSkewSymmetricMatrix(Eigen::Vector3d(
    added_mass_coeff_(0) * velocity(0), added_mass_coeff_(1) * velocity(1), added_mass_coeff_(2) * velocity(2)));

  Eigen::Matrix3d angular_vel = makeSkewSymmetricMatrix(Eigen::Vector3d(
    added_mass_coeff_(3) * velocity(3), added_mass_coeff_(4) * velocity(4), added_mass_coeff_(5) * velocity(5)));

  coriolis.topRightCorner(3, 3) = linear_vel;
  coriolis.bottomLeftCorner(3, 3) = linear_vel;
  coriolis.bottomRightCorner(3, 3) = angular_vel;

  return coriolis;
}

Eigen::Matrix6d Coriolis::calculateCoriolisMatrix(const Eigen::Vector6d & velocity) const
{
  return calculateRigidBodyCoriolisMatrix(velocity.bottomRows(3)) + calculateAddedCoriolisMatrix(velocity);
}

Damping::Damping(double Xu, double Yv, double Zw, double Kp, double Mq, double Nr, double Xuu, double Yvv, double Zww,
                 double Kpp, double Mqq, double Nrr)
: linear_mat_(Eigen::Vector6d(Xu, Yv, Zw, Kp, Mq, Nr).asDiagonal().toDenseMatrix()),
  quadratic_mat_(Eigen::Vector6d(Xuu, Yvv, Zww, Kpp, Mqq, Nrr).asDiagonal().toDenseMatrix())
{
}

Damping::Damping(Eigen::Vector6d linear, Eigen::Vector6d quadratic)
: linear_mat_(linear.asDiagonal().toDenseMatrix()), quadratic_mat_(quadratic.asDiagonal().toDenseMatrix())
{
}

Damping::Damping(Eigen::Matrix6d linear, Eigen::Matrix6d quadratic)
: linear_mat_(std::move(linear)), quadratic_mat_(std::move(quadratic))
{
}

Eigen::Matrix6d Damping::calculateDampingMatrix(const Eigen::Vector6d & velocity) const
{
  Eigen::Matrix6d quadratic = -(quadratic_mat_ * velocity.cwiseAbs()).asDiagonal().toDenseMatrix();
  return -linear_mat_ + quadratic;
}

RestoringForces::RestoringForces(double weight, double buoyancy, Eigen::Vector3d center_of_buoyancy,
                                 Eigen::Vector3d center_of_gravity)
: weight_(weight),
  buoyancy_(buoyancy),
  center_of_buoyancy_(std::move(center_of_buoyancy)),
  center_of_gravity_(std::move(center_of_gravity))
{
}

Eigen::Vector6d RestoringForces::calculateRestoringForcesVector(const Eigen::Matrix3d & rot) const
{
  const Eigen::Vector3d fg(0, 0, weight_);
  const Eigen::Vector3d fb(0, 0, -buoyancy_);

  Eigen::Vector6d g_rb;

  g_rb.topRows(3) = rot * (fg + fb);
  g_rb.bottomRows(3) = center_of_gravity_.cross(rot * fg) + center_of_buoyancy_.cross(rot * fb);

  g_rb *= -1;

  return g_rb;
}

LinkDynamics::LinkDynamics(Inertia inertia, Coriolis coriolis, Damping damping, RestoringForces restoring_forces)
: inertia(std::move(inertia)),
  coriolis(std::move(coriolis)),
  damping(std::move(damping)),
  restoring_forces(std::move(restoring_forces))
{
}

}  // namespace hydrodynamics
