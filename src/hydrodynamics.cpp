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

#include <tinyxml2.h>

#include <Eigen/Dense>
#include <expected>
#include <format>
#include <fstream>
#include <optional>
#include <ranges>
#include <sstream>

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
  const Eigen::Matrix3d & I,
  const Eigen::Matrix6d & added_mass,
  const Eigen::Vector3d & cog)
{
  rigid_body_matrix.topLeftCorner(3, 3) = mass * Eigen::Matrix3d::Identity();
  rigid_body_matrix.topRightCorner(3, 3) = -mass * make_skew_symmetric_matrix(cog);
  rigid_body_matrix.bottomLeftCorner(3, 3) = mass * make_skew_symmetric_matrix(cog);
  rigid_body_matrix.bottomRightCorner(3, 3) = I;
  added_mass_matrix = -added_mass;
  mass_matrix = rigid_body_matrix + added_mass_matrix;
}

Coriolis::Coriolis(
  double mass,
  const Eigen::Matrix3d & I,
  const Eigen::Matrix6d & added_mass,
  const Eigen::Vector3d & cog)
: inertia(mass, I, added_mass, cog)
{
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
auto Coriolis::rigid_body_coriolis_matrix(const Eigen::Vector6d & velocity) const -> Eigen::Matrix6d
{
  const Eigen::Vector3d v1 = velocity.head(3);
  const Eigen::Vector3d v2 = velocity.tail(3);

  const Eigen::Matrix3d m11 = inertia.rigid_body_matrix.topLeftCorner(3, 3);
  const Eigen::Matrix3d m12 = inertia.rigid_body_matrix.topRightCorner(3, 3);
  const Eigen::Matrix3d m21 = inertia.rigid_body_matrix.bottomLeftCorner(3, 3);
  const Eigen::Matrix3d m22 = inertia.rigid_body_matrix.bottomRightCorner(3, 3);

  const Eigen::Vector3d c1 = (m11 * v1) + (m12 * v2);
  const Eigen::Vector3d c2 = (m21 * v1) + (m22 * v2);

  const Eigen::Matrix3d c1_skew = make_skew_symmetric_matrix(c1);
  const Eigen::Matrix3d c2_skew = make_skew_symmetric_matrix(c2);

  Eigen::Matrix6d coriolis = Eigen::Matrix6d::Zero();
  coriolis.topRightCorner(3, 3) = -c1_skew;
  coriolis.bottomLeftCorner(3, 3) = -c1_skew;
  coriolis.bottomRightCorner(3, 3) = -c2_skew;

  return coriolis;
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
auto Coriolis::added_coriolis_matrix(const Eigen::Vector6d & velocity) const -> Eigen::Matrix6d
{
  const Eigen::Vector3d v1 = velocity.head(3);
  const Eigen::Vector3d v2 = velocity.tail(3);

  const Eigen::Matrix3d a11 = inertia.added_mass_matrix.topLeftCorner(3, 3);
  const Eigen::Matrix3d a12 = inertia.added_mass_matrix.topRightCorner(3, 3);
  const Eigen::Matrix3d a21 = inertia.added_mass_matrix.bottomLeftCorner(3, 3);
  const Eigen::Matrix3d a22 = inertia.added_mass_matrix.bottomRightCorner(3, 3);

  const Eigen::Vector3d c1 = (a11 * v1) + (a12 * v2);
  const Eigen::Vector3d c2 = (a21 * v1) + (a22 * v2);

  const Eigen::Matrix3d c1_skew = make_skew_symmetric_matrix(c1);
  const Eigen::Matrix3d c2_skew = make_skew_symmetric_matrix(c2);

  Eigen::Matrix6d coriolis = Eigen::Matrix6d::Zero();
  coriolis.topRightCorner(3, 3) = -c1_skew;
  coriolis.bottomLeftCorner(3, 3) = -c1_skew;
  coriolis.bottomRightCorner(3, 3) = -c2_skew;

  return coriolis;
}

auto Coriolis::coriolis_matrix(const Eigen::Vector6d & velocity) const -> Eigen::Matrix6d
{
  return rigid_body_coriolis_matrix(velocity) + added_coriolis_matrix(velocity);
}

Damping::Damping(const Eigen::Matrix6d & linear, const Eigen::Matrix6d & quadratic)
: linear_drag(-linear),
  quadratic_drag(-quadratic)
{
}

Damping::Damping(Eigen::Vector6d linear, Eigen::Vector6d quadratic)
: linear_drag(linear.asDiagonal().toDenseMatrix()),
  quadratic_drag(quadratic.asDiagonal().toDenseMatrix())
{
}

auto Damping::damping_matrix(const Eigen::Vector6d & velocity) const -> Eigen::Matrix6d
{
  return linear_drag + (quadratic_drag * velocity.cwiseAbs().asDiagonal().toDenseMatrix());
}

RestoringForces::RestoringForces(double weight, double buoyancy, Eigen::Vector3d cob, Eigen::Vector3d cog)
: weight(weight),
  buoyancy(buoyancy),
  center_of_buoyancy(std::move(cob)),
  center_of_gravity(std::move(cog))
{
}

auto RestoringForces::restoring_forces_vector(const Eigen::Matrix3d & rotation) const -> Eigen::Vector6d
{
  const Eigen::Vector3d fg(0, 0, weight);
  const Eigen::Vector3d fb(0, 0, -buoyancy);

  const Eigen::Matrix3d rot_t = rotation.transpose();

  Eigen::Vector6d g_rb;
  g_rb.topRows(3) = rot_t * (fg + fb);
  g_rb.bottomRows(3) = center_of_gravity.cross(rot_t * fg) + center_of_buoyancy.cross(rot_t * fb);
  g_rb *= -1;

  return g_rb;
}

auto RestoringForces::restoring_forces_vector(const Eigen::Quaterniond & rotation) const -> Eigen::Vector6d
{
  return restoring_forces_vector(rotation.toRotationMatrix());
}

Params::Params(const Inertia & M, const Coriolis & C, const Damping & D, const RestoringForces & g)
{
  this->M = M;
  this->C = C;
  this->D = D;
  this->g = g;
}

Params::Params(Inertia && M, Coriolis && C, Damping && D, RestoringForces && g)
{
  this->M = std::move(M);
  this->C = std::move(C);
  this->D = std::move(D);
  this->g = std::move(g);
}

Params::Params(
  double mass,
  const Eigen::Matrix3d & I,
  const Eigen::Matrix6d & added_mass,
  const Eigen::Matrix6d & linear_damping,
  const Eigen::Matrix6d & quadratic_damping,
  const Eigen::Vector3d & cog,
  const Eigen::Vector3d & cob,
  double weight,
  double buoyancy)
: M(mass, I, added_mass, cog),
  C(mass, I, added_mass, cog),
  D(linear_damping, quadratic_damping),
  g(weight, buoyancy, cob, cog)
{
}

// NOLINTNEXTLINE(misc-use-internal-linkage)
auto forward_dynamics(
  const Params & params,
  const Eigen::Vector6d & vel,
  const Eigen::Vector6d & tau,
  const Eigen::Matrix3d & R) -> Eigen::Vector6d
{
  return params.M.mass_matrix.inverse() * (tau - (params.C(vel) * vel) - (params.D(vel) * vel) - params.g(R));
}

// NOLINTNEXTLINE(misc-use-internal-linkage)
auto inverse_dynamics(
  const Params & params,
  const Eigen::Vector6d & acc,
  const Eigen::Vector6d & vel,
  const Eigen::Matrix3d & R) -> Eigen::Vector6d
{
  return (params.M * acc) + (params.C(vel) * vel) + (params.D(vel) * vel) + params.g(R);
}

// anonymous namespace for parsing functions
namespace
{

auto parse_element(const tinyxml2::XMLElement * parent, const std::string & name)
  -> std::expected<const tinyxml2::XMLElement *, std::string>
{
  const tinyxml2::XMLElement * element = parent->FirstChildElement(name.c_str());
  if (element == nullptr) {
    return std::unexpected(std::format("Failed to parse required element: {}", name));
  }
  return element;
}

auto parse_text(const tinyxml2::XMLElement * parent, const std::string & name) -> std::expected<double, std::string>
{
  const auto out = parse_element(parent, name);
  if (!out.has_value()) {
    return std::unexpected(out.error());
  }
  const tinyxml2::XMLElement * element = out.value();
  double value = 0.0;
  if (element->QueryDoubleText(&value) != tinyxml2::XML_SUCCESS) {
    return std::unexpected(std::format("Failed to parse required element: <{}>", name));
  }
  return value;
}

auto parse_attribute(const tinyxml2::XMLElement * element, const std::string & attribute_name)
  -> std::expected<double, std::string>
{
  double value = 0.0;
  if (element == nullptr || element->QueryDoubleAttribute(attribute_name.c_str(), &value) != tinyxml2::XML_SUCCESS) {
    return std::unexpected(std::format("Failed to parse required attribute: {}", attribute_name));
  }
  return value;
}

auto parse_vector6d(const tinyxml2::XMLElement * parent, const std::string & child, std::array<std::string, 6> keys)
  -> std::expected<Eigen::Vector6d, std::string>
{
  Eigen::Vector6d vec = Eigen::Vector6d::Zero();
  const auto out = parse_element(parent, child);
  if (!out.has_value()) {
    return std::unexpected(out.error());
  }
  const tinyxml2::XMLElement * vector_it = out.value();
  for (const auto [i, key] : std::views::enumerate(keys)) {
    const auto val = parse_attribute(vector_it, key);
    if (!val.has_value()) {
      return std::unexpected(val.error());
    }
    vec(i) = val.value();
  }
  return vec;
}

auto parse_vector3d(const tinyxml2::XMLElement * parent, const std::string & child)
  -> std::expected<Eigen::Vector3d, std::string>
{
  Eigen::Vector3d vec = Eigen::Vector3d::Zero();
  const auto out = parse_element(parent, child);
  if (!out.has_value()) {
    return std::unexpected(out.error());
  }
  const tinyxml2::XMLElement * vector_it = out.value();
  vector_it->QueryDoubleAttribute("x", &vec(0));
  vector_it->QueryDoubleAttribute("y", &vec(1));
  vector_it->QueryDoubleAttribute("z", &vec(2));
  return vec;
}

}  // namespace

// NOLINTNEXTLINE(misc-use-internal-linkage)
[[nodiscard]] auto parse_model_from_xml(const std::string & tree) -> std::expected<Params, std::string>
{
  if (tree.empty()) {
    return std::unexpected("URDF tree is empty");
  }

  tinyxml2::XMLDocument doc;
  if (!doc.Parse(tree.c_str()) && doc.Error()) {
    return std::unexpected("Failed to parse URDF tree");
  }

  const tinyxml2::XMLElement * robot_it = doc.RootElement();
  if (robot_it == nullptr || robot_it->Name() != std::string("robot")) {
    return std::unexpected("Invalid URDF tree: missing <robot> root");
  }

  const tinyxml2::XMLElement * hydro_it = robot_it->FirstChildElement("hydrodynamics");
  if (hydro_it == nullptr) {
    return std::unexpected("No <hydrodynamics> element found in URDF tree");
  }

  // extract the mass, weight, and buoyancy of the vehicle
  std::unordered_map<std::string, double> text_params;
  for (const auto key : {"mass", "weight", "buoyancy"}) {
    const auto val = parse_text(hydro_it, key);
    if (!val.has_value()) {
      return std::unexpected(val.error());
    }
    text_params[key] = val.value();
  }
  // NOLINTNEXTLINE
  const double mass = text_params["mass"], weight = text_params["weight"], buoyancy = text_params["buoyancy"];

  // extract the inertia tensor
  const auto inertia_out = parse_element(hydro_it, "inertia");
  if (!inertia_out.has_value()) {
    return std::unexpected(inertia_out.error());
  }
  const tinyxml2::XMLElement * inertia_it = inertia_out.value();

  // the moments of inertia are required
  Eigen::Matrix3d inertia = Eigen::Matrix3d::Zero();
  const std::array<std::string, 3> inertia_keys = {"ixx", "iyy", "izz"};
  for (const auto [i, key] : std::views::enumerate(inertia_keys)) {
    const auto val = parse_attribute(inertia_it, key);
    if (!val.has_value()) {
      return std::unexpected(val.error());
    }
    inertia(i, i) = val.value();
  }

  // the products of inertia are optional
  inertia(0, 1) = parse_attribute(inertia_it, "ixy").value_or(0.0);
  inertia(0, 2) = parse_attribute(inertia_it, "ixz").value_or(0.0);
  inertia(1, 2) = parse_attribute(inertia_it, "iyz").value_or(0.0);

  // parse the added mass coefficients
  const auto added_mass_out = parse_vector6d(hydro_it, "added_mass", {"Xdu", "Ydv", "Zdw", "Kdp", "Mdq", "Ndr"});
  if (!added_mass_out.has_value()) {
    return std::unexpected(added_mass_out.error());
  }
  Eigen::Matrix6d added_mass = added_mass_out.value().asDiagonal().toDenseMatrix();

  // parse the linear and quadratic damping coefficients
  const auto linear_drag_out = parse_vector6d(hydro_it, "linear_damping", {"Xu", "Yv", "Zw", "Kp", "Mq", "Nr"});
  if (!linear_drag_out.has_value()) {
    return std::unexpected(linear_drag_out.error());
  }
  Eigen::Matrix6d linear_drag = linear_drag_out.value().asDiagonal().toDenseMatrix();

  const auto quad_drag_out = parse_vector6d(hydro_it, "quadratic_damping", {"Xuu", "Yvv", "Zww", "Kpp", "Mqq", "Nrr"});
  if (!quad_drag_out.has_value()) {
    return std::unexpected(quad_drag_out.error());
  }
  Eigen::Matrix6d quadratic_drag = quad_drag_out.value().asDiagonal().toDenseMatrix();

  // parse the center of gravity
  Eigen::Vector3d cog = Eigen::Vector3d::Zero();
  const auto cog_out = parse_vector3d(hydro_it, "center_of_gravity");
  if (!cog_out.has_value()) {
    return std::unexpected(cog_out.error());
  }
  cog = cog_out.value();

  // parse the center of buoyancy
  Eigen::Vector3d cob = Eigen::Vector3d::Zero();
  const auto cob_out = parse_vector3d(hydro_it, "center_of_buoyancy");
  if (!cob_out.has_value()) {
    return std::unexpected(cob_out.error());
  }
  cob = cob_out.value();

  Params params(mass, inertia, added_mass, linear_drag, quadratic_drag, cog, cob, weight, buoyancy);
  return params;
}

// NOLINTNEXTLINE(misc-use-internal-linkage)
[[nodiscard]] auto parse_model_from_urdf(const std::string & file) -> std::expected<Params, std::string>
{
  if (file.empty()) {
    return std::unexpected("URDF file path is empty");
  }
  const std::ifstream buffer(file);
  std::stringstream ss;
  ss << buffer.rdbuf();
  return parse_model_from_xml(ss.str());
}

}  // namespace hydrodynamics
