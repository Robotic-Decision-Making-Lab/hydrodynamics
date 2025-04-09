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

#include <format>
#include <iostream>
#include <string>

#include "hydrodynamics/hydrodynamics.hpp"

auto main() -> int
{
  // parse the urdf from a file path - update the path to the urdf file
  // we use the BlueROV2 hydrodynamics model as an example
  const std::string urdf_path = "/home/ubuntu/ws_ros/src/hydrodynamics/examples/description/bluerov2.model.urdf";
  const auto out = hydrodynamics::parse_model_from_urdf(urdf_path);

  // check if the parsing was successful, if not print the error message
  if (!out.has_value()) {
    std::cerr << std::format("Failed to parse URDF: {}\n", out.error());
    return 1;
  }
  std::cout << "Parsed URDF successfully!\n";

  // get the parsed parameters
  const hydrodynamics::Params model = out.value();

  const Eigen::Vector6d acc = Eigen::Vector6d::Zero();
  const Eigen::Vector6d vel = Eigen::Vector6d::Ones();
  const Eigen::Matrix3d rot = Eigen::Matrix3d::Identity();

  // perform inverse dynamics using the parsed parameters
  const Eigen::Vector6d tau = hydrodynamics::inverse_dynamics(model, acc, vel, rot);

  std::cout << "Result of inverse dynamics:\n";
  std::cout << std::format(
    "fx: {}, fy: {}, fz: {}, tx: {}, ty: {}, tz: {}\n", tau(0), tau(1), tau(2), tau(3), tau(4), tau(5));

  return 0;
}
