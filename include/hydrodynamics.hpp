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

#include "eigen.hpp"

namespace hydrodynamics
{

/**
 * @brief A wrapper for inertial parameters including added mass and rigid body inertia.
 */
class Inertia
{
public:
  /**
   * @brief Construct a new wrapper for inertial parameters including added mass and rigid body inertia.
   *
   * @param mass The mass of the vehicle.
   * @param Ixx The moment of inertia about the x-axis.
   * @param Iyy The moment of inertia about the y-axis.
   * @param Izz The moment of inertia about the z-axis.
   * @param Xdu The added mass coefficient in the surge direction.
   * @param Ydv The added mass coefficient in the sway direction.
   * @param Zdw The added mass coefficient in the heave direction.
   * @param Kdp The added mass coefficient in the roll direction.
   * @param Mdq The added mass coefficient in the pitch direction.
   * @param Ndr The added mass coefficient in the yaw direction.
   */
  Inertia(double mass, double Ixx, double Iyy, double Izz, double Xdu, double Ydv, double Zdw, double Kdp, double Mdq,
          double Ndr);

  /**
   * @brief Construct a new wrapper for inertial parameters including added mass and rigid body inertia.
   *
   * @param mass The mass of the vehicle.
   * @param moments The moments of inertia about the x, y, and z axes.
   * @param added_mass The added mass coefficients in the surge, sway, heave, roll, pitch, and yaw directions.
   */
  Inertia(double mass, Eigen::Vector3d moments, Eigen::Vector6d added_mass);

  /**
   * @brief Get the rigid body inertia matrix.
   *
   * @return Eigen::Matrix6d
   */
  [[nodiscard]] Eigen::Matrix6d getRigidBodyInertiaMatrix() const;

  /**
   * @brief Get the added mass matrix.
   *
   * @return Eigen::Matrix6d
   */
  [[nodiscard]] Eigen::Matrix6d getAddedMassMatrix() const;

  /**
   * @brief Get the mass matrix.
   *
   * @note This is the sum of the rigid body and added mass matrices.
   *
   * @return Eigen::Matrix6d
   */
  [[nodiscard]] Eigen::Matrix6d getMassMatrix() const;

private:
  Eigen::Matrix6d rigid_body_mat_;
  Eigen::Matrix6d added_mass_mat_;
  Eigen::Matrix6d mass_mat_;
};

/**
 * @brief A wrapper for Coriolis and centripetal force parameters.
 */
class Coriolis
{
public:
  /**
   * @brief Construct a new wrapper for the Coriolis and centripetal force parameters.
   *
   * @param mass The mass of the vehicle.
   * @param Ixx The moment of inertia about the x-axis.
   * @param Iyy The moment of inertia about the y-axis.
   * @param Izz The moment of inertia about the z-axis.
   * @param Xdu The added mass coefficient in the surge direction.
   * @param Ydv The added mass coefficient in the sway direction.
   * @param Zdw The added mass coefficient in the heave direction.
   * @param Kdp The added mass coefficient in the roll direction.
   * @param Mdq The added mass coefficient in the pitch direction.
   * @param Ndr The added mass coefficient in the yaw direction.
   */
  Coriolis(double mass, double Ixx, double Iyy, double Izz, double Xdu, double Ydv, double Zdw, double Kdp, double Mdq,
           double Ndr);

  /**
   * @brief Construct a new wrapper for the Coriolis and centripetal force parameters.
   *
   * @param mass The mass of the vehicle.
   * @param moments The moments of inertia about the x, y, and z axes.
   * @param added_mass The added mass coefficients in the surge, sway, heave, roll, pitch, and yaw directions.
   */
  Coriolis(double mass, Eigen::Vector3d moments, Eigen::Vector6d added_mass);

  /**
   * @brief Calculate the rigid body Coriolis matrix.
   *
   * @param angular_velocity The angular velocity of the vehicle.
   * @return Eigen::Matrix6d
   */
  [[nodiscard]] Eigen::Matrix6d calculateRigidBodyCoriolisMatrix(const Eigen::Vector3d & angular_velocity) const;

  /**
   * @brief Calculate the added mass Coriolis matrix.
   *
   * @param velocity The velocity of the vehicle.
   * @return Eigen::Matrix6d
   */
  [[nodiscard]] Eigen::Matrix6d calculateAddedCoriolisMatrix(const Eigen::Vector6d & velocity) const;

  /**
   * @brief Calculate the Coriolis matrix.
   *
   * @note This is the sum of the rigid body and added mass Coriolis matrices.
   *
   * @param velocity The velocity of the vehicle.
   * @return Eigen::Matrix6d
   */
  [[nodiscard]] Eigen::Matrix6d calculateCoriolisMatrix(const Eigen::Vector6d & velocity) const;

private:
  double mass_;
  Eigen::Matrix3d moments_;
  Eigen::Vector6d added_mass_coeff_;
};

/**
 * @brief A wrapper for the damping forces acting on the vehicle, including linear and quadratic damping coefficients.
 */
class Damping
{
public:
  /**
   * @brief Construct a new wrapper for the damping coefficients.
   *
   * @param Xu The linear damping coefficient in the surge direction.
   * @param Yv The linear damping coefficient in the sway direction.
   * @param Zw The linear damping coefficient in the heave direction.
   * @param Kp The linear damping coefficient in the roll direction.
   * @param Mq The linear damping coefficient in the pitch direction.
   * @param Nr The linear damping coefficient in the yaw direction.
   * @param Xuu The quadratic damping coefficient in the surge direction.
   * @param Yvv The quadratic damping coefficient in the sway direction.
   * @param Zww The quadratic damping coefficient in the heave direction.
   * @param Kpp The quadratic damping coefficient in the roll direction.
   * @param Mqq The quadratic damping coefficient in the pitch direction.
   * @param Nrr The quadratic damping coefficient in the yaw direction.
   */
  Damping(double Xu, double Yv, double Zw, double Kp, double Mq, double Nr, double Xuu, double Yvv, double Zww,
          double Kpp, double Mqq, double Nrr);

  /**
   * @brief Construct a new wrapper for the damping coefficients.
   *
   * @param linear The linear damping coefficients.
   * @param quadratic The quadratic damping coefficients.
   */
  Damping(Eigen::Vector6d linear, Eigen::Vector6d quadratic);

  /**
   * @brief Calculate the damping matrix.
   *
   * @param velocity The current velocity of the vehicle.
   * @return Eigen::Matrix6d
   */
  [[nodiscard]] Eigen::Matrix6d calculateDampingMatrix(const Eigen::Vector6d & velocity) const;

private:
  Eigen::Vector6d linear_coeff_;
  Eigen::Vector6d quadratic_coeff_;
};

/**
 * @brief A wrapper for the restoring forces acting on the vehicle, including gravity and bouyancy.
 */
class RestoringForces
{
public:
  /**
   * @brief Construct a wrapper for restoring force parameters.
   *
   * @param weight The weight of the vehicle.
   * @param buoyancy The buoyancy of the vehicle.
   * @param center_of_buoyancy The center of buoyancy of the vehicle.
   * @param center_of_gravity The center of gravity of the vehicle.
   */
  RestoringForces(double weight, double buoyancy, Eigen::Vector3d center_of_buoyancy,
                  Eigen::Vector3d center_of_gravity);

  /**
   * @brief Calculate the restoring forces acting on the vehicle.
   *
   * @param rotation The current rotation of the vehicle with respect to the inertial frame.
   * @return Eigen::Vector6d
   */
  [[nodiscard]] Eigen::Vector6d calculateRestoringForcesVector(const Eigen::Matrix3d & rotation) const;

private:
  double weight_;
  double buoyancy_;
  Eigen::Vector3d center_of_buoyancy_;
  Eigen::Vector3d center_of_gravity_;
};

}  // namespace hydrodynamics
