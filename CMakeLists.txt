cmake_minimum_required(VERSION 3.8)
project(hydrodynamics)

if(CMAKE_CXX_COMPILER_ID MATCHES "(GNU|Clang)")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

set(THIS_PACKAGE_INCLUDE_DEPENDS
  rclcpp
  Eigen3
)

find_package(ament_cmake REQUIRED)
foreach(Dependency IN ITEMS ${THIS_PACKAGE_INCLUDE_DEPENDS})
  find_package(${Dependency} REQUIRED)
endforeach()

add_library(hydrodynamics SHARED
  src/hydrodynamics.cpp
)

target_include_directories(hydrodynamics
  PUBLIC
    $<INSTALL_INTERFACE:include/hydrodynamics>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    ${EIGEN3_INCLUDE_DIR}
)
target_compile_features(hydrodynamics PUBLIC cxx_std_17)
target_link_libraries(hydrodynamics
  PUBLIC
    ${rclcpp_LIBRARIES}
)

ament_target_dependencies(hydrodynamics
  PUBLIC
    rclcpp
)

install(
  DIRECTORY include/
  DESTINATION include/hydrodynamics
)

install(
  TARGETS hydrodynamics
  EXPORT export_hydrodynamics
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION include
)

if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)

  # Run linters found in package.xml except those below
  set(ament_cmake_copyright_FOUND TRUE)
  set(ament_cmake_uncrustify_FOUND TRUE)
  set(ament_cmake_pep257_FOUND TRUE)
  set(ament_cmake_flake8_FOUND TRUE)
endif()

ament_export_targets(export_hydrodynamics HAS_LIBRARY_TARGET)
ament_export_dependencies(${THIS_PACKAGE_INCLUDE_DEPENDS})

ament_package()
