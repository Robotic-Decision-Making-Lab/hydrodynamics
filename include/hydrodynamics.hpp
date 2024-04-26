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
 * @brief This class provides functionality to store and compute the system inertia matrix for a submerged body.
 *
 * @note The inertia matrix is based off of Fossen's equations of motion for marine vehicles [1] and includes both the
 * rigid body inertia matrix and the added mass matrix.
 *
 * [1] Fossen, Thor I. _Guidance and Control of Ocean Vehicles_.
 *     United Kingdom: Wiley, 1994.
 */
class Inertia
{
public:
  /**
   * @brief Constructor.
   *
   * @param mass Mass of the body [kg].
   * @param Ixx Moment about the x-axis [kgm^2].
   * @param Iyy Moment about the y-axis [kgm^2].
   * @param Izz Moment about the z-axis [kgm^2].
   * @param Xdu Added mass in the surge direction [kg].
   * @param Ydv Added mass in the sway direction [kg].
   * @param Zdw Added mass in the heave direction [kg].
   * @param Kdp Added mass in the roll direction [kgm^2].
   * @param Mdq Added mass in the pitch direction [kgm^2].
   * @param Ndr Added mass in the yaw direction [kgm^2].
   */
  Inertia(double mass, double Ixx, double Iyy, double Izz, double Xdu, double Ydv, double Zdw, double Kdp, double Mdq,
          double Ndr);

  /**
   * @brief Constructor.
   *
   * @param mass Mass of the body [kg].
   * @param Ixx Moment about the x-axis [kgm^2].
   * @param Ixy Coupled moment about the x and y axes [kgm^2].
   * @param Ixz Coupled moment about the x and z axes [kgm^2].
   * @param Iyy Moment about the y-axis [kgm^2].
   * @param Iyz Coupled moment about the y and z axes [kgm^2].
   * @param Izz Moment about the z-axis [kgm^2].
   * @param Xdu Added mass in the surge direction [kg].
   * @param Ydv Added mass in the sway direction [kg].
   * @param Zdw Added mass in the heave direction [kg].
   * @param Kdp Added mass in the roll direction [kgm^2].
   * @param Mdq Added mass in the pitch direction [kgm^2].
   * @param Ndr Added mass in the yaw direction [kgm^2].
   */
  Inertia(double mass, double Ixx, double Ixy, double Ixz, double Iyy, double Iyz, double Izz, double Xdu, double Ydv,
          double Zdw, double Kdp, double Mdq, double Ndr);

  /**
   * @brief Constructor.
   *
   * @param mass Mass of the body [kg].
   * @param moments Moments about the x, y, and z axes, respectively [kgm^2].
   * @param added_mass Added mass coefficients in the surge, sway, heave, roll, pitch, and yaw directions [kg, kgm^2].
   */
  Inertia(double mass, Eigen::Vector3d moments, Eigen::Vector6d added_mass);

  /**
   * @brief Constructor.
   *
   * @param mass Mass of the body [kg].
   * @param moments Moments of inertia [kgm^2].
   * @param added_mass Added mass coefficients in the surge, sway, heave, roll, pitch, and yaw directions [kg, kgm^2].
   */
  Inertia(double mass, Eigen::Matrix3d moments, Eigen::Vector6d added_mass);

  /**
   * @brief Get the rigid body inertia matrix.
   *
   * @return Eigen::Matrix6d
   */
  [[nodiscard]] Eigen::Matrix6d getRigidBodyInertiaMatrix() const { return rigid_body_mat_; }

  /**
   * @brief Get the added mass matrix.
   *
   * @return Eigen::Matrix6d
   */
  [[nodiscard]] Eigen::Matrix6d getAddedMassMatrix() const { return added_mass_mat_; };

  /**
   * @brief Get the mass matrix.
   *
   * @note This is the sum of the rigid body and added mass matrices.
   *
   * @return Eigen::Matrix6d
   */
  [[nodiscard]] Eigen::Matrix6d getSystemInertiaMatrix() const { return inertia_mat_; }

private:
  Eigen::Matrix6d rigid_body_mat_;
  Eigen::Matrix6d added_mass_mat_;
  Eigen::Matrix6d inertia_mat_;
};

/**
 * @brief This class provides functionality to store and compute the system Coriolis matrix for a submerged body.
 *
 * @note The Coriolis matrix is based off of Fossen's equations of motion for marine vehicles [1] and includes both the
 * rigid body Coriolis matrix and the added mass Coriolis matrix.
 *
 * [1] Fossen, Thor I. _Guidance and Control of Ocean Vehicles_.
 *     United Kingdom: Wiley, 1994.
 */
class Coriolis
{
public:
  /**
   * @brief Constructor.
   *
   * @param mass Mass of the body [kg].
   * @param Ixx Moment about the x-axis [kgm^2].
   * @param Iyy Moment about the y-axis [kgm^2].
   * @param Izz Moment about the z-axis [kgm^2].
   * @param Xdu Added mass in the surge direction [kg].
   * @param Ydv Added mass in the sway direction [kg].
   * @param Zdw Added mass in the heave direction [kg].
   * @param Kdp Added mass in the roll direction [kgm^2].
   * @param Mdq Added mass in the pitch direction [kgm^2].
   * @param Ndr Added mass in the yaw direction [kgm^2].
   */
  Coriolis(double mass, double Ixx, double Iyy, double Izz, double Xdu, double Ydv, double Zdw, double Kdp, double Mdq,
           double Ndr);

  /**
   * @brief Constructor.
   *
   * @param mass Mass of the body [kg].
   * @param Ixx Moment about the x-axis [kgm^2].
   * @param Ixy Coupled moment about the x and y axes [kgm^2].
   * @param Ixz Coupled moment about the x and z axes [kgm^2].
   * @param Iyy Moment about the y-axis [kgm^2].
   * @param Iyz Coupled moment about the y and z axes [kgm^2].
   * @param Izz Moment about the z-axis [kgm^2].
   * @param Xdu Added mass in the surge direction [kg].
   * @param Ydv Added mass in the sway direction [kg].
   * @param Zdw Added mass in the heave direction [kg].
   * @param Kdp Added mass in the roll direction [kgm^2].
   * @param Mdq Added mass in the pitch direction [kgm^2].
   * @param Ndr Added mass in the yaw direction [kgm^2].
   */
  Coriolis(double mass, double Ixx, double Ixy, double Ixz, double Iyy, double Iyz, double Izz, double Xdu, double Ydv,
           double Zdw, double Kdp, double Mdq, double Ndr);

  /**
   * @brief Constructor.
   *
   * @param mass Mass of the body [kg].
   * @param moments Moments about the x, y, and z axes, respectively [kgm^2].
   * @param added_mass Added mass coefficients in the surge, sway, heave, roll, pitch, and yaw directions [kg, kgm^2].
   */
  Coriolis(double mass, Eigen::Vector3d moments, Eigen::Vector6d added_mass);

  /**
   * @brief Constructor.
   *
   * @param mass Mass of the body [kg].
   * @param moments Moments of inertia [kgm^2].
   * @param added_mass Added mass coefficients in the surge, sway, heave, roll, pitch, and yaw directions [kg, kgm^2].
   */
  Coriolis(double mass, Eigen::Matrix3d moments, Eigen::Vector6d added_mass);

  /**
   * @brief Calculate the rigid body Coriolis matrix.
   *
   * @param angular_velocity The angular velocity of the body.
   * @return Eigen::Matrix6d
   */
  [[nodiscard]] Eigen::Matrix6d calculateRigidBodyCoriolisMatrix(const Eigen::Vector3d & angular_velocity) const;

  /**
   * @brief Calculate the added mass Coriolis matrix.
   *
   * @param velocity The velocity of the body.
   * @return Eigen::Matrix6d
   */
  [[nodiscard]] Eigen::Matrix6d calculateAddedCoriolisMatrix(const Eigen::Vector6d & velocity) const;

  /**
   * @brief Calculate the Coriolis matrix.
   *
   * @note This is the sum of the rigid body and added mass Coriolis matrices.
   *
   * @param velocity The velocity of the body.
   * @return Eigen::Matrix6d
   */
  [[nodiscard]] Eigen::Matrix6d calculateCoriolisMatrix(const Eigen::Vector6d & velocity) const;

private:
  double mass_;
  Eigen::Matrix3d moments_mat_;
  Eigen::Vector6d added_mass_coeff_;
};

/**
 * @brief This class provides functionality to store and compute the viscous damping matrix for a submerged body.
 *
 * @note The damping matrix is based off of Fossen's equations of motion for marine vehicles [1] and includes
 * linear and quadratic damping.
 *
 * [1] Fossen, Thor I. _Guidance and Control of Ocean Vehicles_.
 *     United Kingdom: Wiley, 1994.
 */
class Damping
{
public:
  /**
   * @brief Constructor.
   *
   * @param Xu Stability derivative, 1st order, surge component [kg]
   * @param Yv Stability derivative, 1st order, sway component [kg]
   * @param Zw Stability derivative, 1st order, heave component [kg]
   * @param Kp Stability derivative, 1st order, roll component [kg/m]
   * @param Mq Stability derivative, 1st order, pitch component [kg/m]
   * @param Nr Stability derivative, 1st order, yaw component [kg/m]
   * @param Xuu Stability derivative, 2nd order, surge component [kg/m]
   * @param Yvv Stability derivative, 2nd order, sway component [kg/m]
   * @param Zww Stability derivative, 2nd order, heave component [kg/m]
   * @param Kpp Stability derivative, 2nd order, roll component [kg/m^2]
   * @param Mqq Stability derivative, 2nd order, pitch component [kg/m^2]
   * @param Nrr Stability derivative, 2nd order, yaw component [kg/m^2]
   */
  Damping(double Xu, double Yv, double Zw, double Kp, double Mq, double Nr, double Xuu, double Yvv, double Zww,
          double Kpp, double Mqq, double Nrr);

  /**
   * @brief Constructor.
   *
   * @param linear Linear damping coefficients in the surge, sway, heave, roll, pitch, and yaw directions [kg, kg/m].
   * @param quadratic Quadratic damping coefficients in the surge, sway, heave, roll, pitch and yaw directions
   *                  [kg/m, kg/m^2].
   */
  Damping(Eigen::Vector6d linear, Eigen::Vector6d quadratic);

  /**
   * @brief Constructor.
   *
   * @param linear Linear damping coefficients, including coupling dissapative terms [kg, kg/m].
   * @param quadratic Quadratic damping coefficients, including coupling dissapative terms [kg/m, kg/m^2].
   */
  Damping(Eigen::Matrix6d linear, Eigen::Matrix6d quadratic);

  /**
   * @brief Calculate the damping matrix.
   *
   * @param velocity The current velocity of the body [m/s].
   * @return Eigen::Matrix6d
   */
  [[nodiscard]] Eigen::Matrix6d calculateDampingMatrix(const Eigen::Vector6d & velocity) const;

private:
  Eigen::Matrix6d linear_mat_;
  Eigen::Matrix6d quadratic_mat_;
};

/**
 * @brief This class provides functionality to store and compute the restoring forces for a submerged body.
 *
 * @note The restoring forces implementation is based off of Fossen's equations of motion for marine vehicles [1], and
 * includes the weight, buoyancy, center of buoyancy, and center of gravity.
 *
 * [1] Fossen, Thor I. _Guidance and Control of Ocean Vehicles_.
 *     United Kingdom: Wiley, 1994.
 */
class RestoringForces
{
public:
  /**
   * @brief Constructor.
   *
   * @param weight Weight of the body (mass * gravity).
   * @param buoyancy Buoyancy of the body (fluid displacement * gravity * water density).
   * @param center_of_buoyancy Center of buoyancy of the body.
   * @param center_of_gravity Center of gravity of the body.
   */
  RestoringForces(double weight, double buoyancy, Eigen::Vector3d center_of_buoyancy,
                  Eigen::Vector3d center_of_gravity);

  /**
   * @brief Calculate the restoring forces acting on the body.
   *
   * @param rotation Current rotation of the body with respect to the inertial frame.
   * @return Eigen::Vector6d
   */
  [[nodiscard]] Eigen::Vector6d calculateRestoringForcesVector(const Eigen::Matrix3d & rotation) const;

private:
  double weight_;
  double buoyancy_;
  Eigen::Vector3d center_of_buoyancy_;
  Eigen::Vector3d center_of_gravity_;
};

/**
 * @brief A wrapper for the dynamic parameters of a submerged body.
 */
struct LinkDynamics
{
  LinkDynamics(Inertia inertia, Coriolis coriolis, Damping damping, RestoringForces restoring_forces);

  Inertia inertia;
  Coriolis coriolis;
  Damping damping;
  RestoringForces restoring_forces;
};

}  // namespace hydrodynamics
