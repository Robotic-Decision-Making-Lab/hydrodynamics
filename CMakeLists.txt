cmake_minimum_required(VERSION 3.8)
project(hydrodynamics)

if(CMAKE_CXX_COMPILER_ID MATCHES "(GNU|Clang)")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

find_package(ament_cmake REQUIRED)
find_package(Eigen3 REQUIRED)

add_library(hydrodynamics SHARED)

target_sources(hydrodynamics
  PRIVATE
    src/hydrodynamics.cpp
)
target_include_directories(hydrodynamics
  PUBLIC
    $<INSTALL_INTERFACE:include/hydrodynamics>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)
target_compile_features(hydrodynamics PUBLIC cxx_std_20)
target_link_libraries(hydrodynamics Eigen3::Eigen)

install(
  DIRECTORY include/
  DESTINATION include/hydrodynamics
)

install(
  TARGETS hydrodynamics
  EXPORT export_hydrodynamics
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
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
ament_export_dependencies(Eigen3)

ament_package()
