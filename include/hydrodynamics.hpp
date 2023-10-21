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

#include "hydrodynamics/eigen.hpp"

namespace hydrodynamics {

struct Inertia {
  Inertia() = default;

  Inertia(double mass, double Ixx, double Iyy, double Izz, double Xdu,
          double Ydv, double Zdw, double Kdp, double Mdq, double Ndr);

  Inertia(double mass, double Ixx, double Iyy, double Izz, double Ixy,
          double Ixz, double Iyz, double Xdu, double Ydv, double Zdw,
          double Kdp, double Mdq, double Ndr);

  Inertia(double mass, const Eigen::Vector3d &moments,
          const Eigen::Vector6d &added_mass);

  Inertia(double mass, const Eigen::Vector6d &moments,
          const Eigen::Vector6d &added_mass);

  Inertia(double mass, const Eigen::Matrix3d &moments,
          const Eigen::Vector6d &added_mass);

  Inertia(double mass, const Eigen::Vector6d &moments,
          const Eigen::Vector6d &added_mass);

  double mass;
  std::shared_ptr<Eigen::Matrix6d> matrix;
  std::shared_ptr<Eigen::Matrix3d> moments;
  std::shared_ptr<Eigen::Matrix6d> added_mass;
};

struct Coriolis {
  Coriolis() = default;

  Coriolis(double mass, double Ixx, double Iyy, double Izz, double Xdu,
           double Ydv, double Zdw, double Kdp, double Mdq, double Ndr);

  Coriolis(double mass, double Ixx, double Iyy, double Izz, double Ixy,
           double Ixz, double Iyz, double Xdu, double Ydv, double Zdw,
           double Kdp, double Mdq, double Ndr);

  Coriolis(double mass, const Eigen::Vector3d &moments,
           const Eigen::Vector6d &added_mass);

  Coriolis(double mass, const Eigen::Vector6d &moments,
           const Eigen::Vector6d &added_mass);

  Coriolis(double mass, const Eigen::Matrix3d &moments,
           const Eigen::Vector6d &added_mass);

  Coriolis(double mass, const Eigen::Vector6d &moments,
           const Eigen::Vector6d &added_mass);

  [[nodiscard]] Eigen::Matrix6d
  calculateRigidBodyCoriolis(const Eigen::Vector3d &angular_velocity) const;

  [[nodiscard]] Eigen::Matrix6d
  calculateAddedMassCoriolis(const Eigen::Vector6d &velocity) const;

  [[nodiscard]] Eigen::Matrix6d
  calculateCoriolis(const Eigen::Vector6d &velocity) const;

  double mass;
  std::shared_ptr<Eigen::Matrix3d> moments;
  std::shared_ptr<Eigen::Matrix6d> added_mass;
};

struct Damping {
  Damping() = default;

  Damping(double Xu, double Yv, double Zw, double Kp, double Mq, double Nr,
          double Xuu, double Yvv, double Zww, double Kpp, double Mqq,
          double Nrr);

  Damping(const Eigen::Vector6d &linear, const Eigen::Vector6d &quadratic);

  [[nodiscard]] Eigen::Matrix6d
  calculateDamping(const Eigen::Vector6d &velocity) const;

  std::shared_ptr<Eigen::Vector6d> linear;
  std::shared_ptr<Eigen::Vector6d> quadratic;
};

struct RestoringForces {
  RestoringForces() = default;

  RestoringForces(double weight, double buoyancy,
                  Eigen::Vector3d center_of_buoyancy,
                  Eigen::Vector3d center_of_gravity);

  [[nodiscard]] Eigen::Vector6d
  calculateRestoringForces(const Eigen::Matrix3d &rot) const;

  double weight;
  double buoyancy;
  std::shared_ptr<Eigen::Vector3d> center_of_buoyancy;
  std::shared_ptr<Eigen::Vector3d> center_of_gravity;
};

struct CurrentEffects {
  CurrentEffects() = default;

  CurrentEffects(Eigen::Vector6d velocity);

  [[nodiscard]] Eigen::Vector6d
  calculateCurrentEffects(const Eigen::Matrix3d &rot) const;

  std::shared_ptr<Eigen::Vector6d> velocity;
};

} // namespace hydrodynamics