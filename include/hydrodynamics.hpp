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

#pragma once

#include <Eigen/Dense>
#include <memory>

#include "eigen.hpp"

namespace hydrodynamics
{

class Inertia
{
public:
  Inertia() = default;

  Inertia(double mass,
          double Ixx,
          double Iyy,
          double Izz,
          double Xdu,
          double Ydv,
          double Zdw,
          double Kdp,
          double Mdq,
          double Ndr);

  Inertia(double mass, const Eigen::Vector3d & moments, const Eigen::Vector6d & added_mass);

  [[nodiscard]] Eigen::Matrix6d getRigidBodyInertiaMatrix() const;

  [[nodiscard]] Eigen::Matrix6d getAddedMassMatrix() const;

  [[nodiscard]] Eigen::Matrix6d getMassMatrix() const;

private:
  Eigen::Matrix6d rigid_body_mat_;
  Eigen::Matrix6d added_mass_mat_;
  Eigen::Matrix6d mass_mat_;
};

class Coriolis
{
public:
  Coriolis() = default;

  Coriolis(double mass,
           double Ixx,
           double Iyy,
           double Izz,
           double Xdu,
           double Ydv,
           double Zdw,
           double Kdp,
           double Mdq,
           double Ndr);

  Coriolis(double mass, const Eigen::Vector3d & moments, Eigen::Vector6d added_mass);

  [[nodiscard]] Eigen::Matrix6d calculateRigidBodyCoriolisMatrix(
    const Eigen::Vector3d & angular_velocity) const;

  [[nodiscard]] Eigen::Matrix6d calculateAddedCoriolisMatrix(
    const Eigen::Vector6d & velocity) const;

  [[nodiscard]] Eigen::Matrix6d calculateCoriolisMatrix(const Eigen::Vector6d & velocity) const;

private:
  double mass_;
  Eigen::Matrix3d moments_;
  Eigen::Vector6d added_mass_coeff_;
};

class Damping
{
public:
  Damping() = default;

  Damping(double Xu,
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

  Damping(Eigen::Vector6d linear, Eigen::Vector6d quadratic);

  [[nodiscard]] Eigen::Matrix6d calculateDampingMatrix(const Eigen::Vector6d & velocity) const;

private:
  Eigen::Vector6d linear_coeff_;
  Eigen::Vector6d quadratic_coeff_;
};

class RestoringForces
{
public:
  RestoringForces() = default;

  RestoringForces(double weight,
                  double buoyancy,
                  Eigen::Vector3d center_of_buoyancy,
                  Eigen::Vector3d center_of_gravity);

  [[nodiscard]] Eigen::Vector6d calculateRestoringForcesVector(const Eigen::Matrix3d & rot) const;

private:
  double weight_;
  double buoyancy_;
  Eigen::Vector3d center_of_buoyancy_;
  Eigen::Vector3d center_of_gravity_;
};

class CurrentEffects
{
  CurrentEffects() = default;

  CurrentEffects(Eigen::Vector6d velocity);

  [[nodiscard]] Eigen::Vector6d calculateCurrentEffectsVector(const Eigen::Matrix3d & rot) const;

private:
  Eigen::Vector6d current_;
};

}  // namespace hydrodynamics
