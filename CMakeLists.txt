cmake_minimum_required(VERSION 3.23)
project(hydrodynamics)

if(CMAKE_CXX_COMPILER_ID MATCHES "(GNU|Clang)")
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

find_package(ament_cmake REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(TinyXML2 REQUIRED)
find_package(tinyxml2_vendor REQUIRED)

add_library(hydrodynamics SHARED)
target_sources(
    hydrodynamics
    PRIVATE src/hydrodynamics.cpp
    PUBLIC
        FILE_SET HEADERS
        BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
        FILES
            ${CMAKE_CURRENT_SOURCE_DIR}/include/hydrodynamics/hydrodynamics.hpp
)

target_compile_features(hydrodynamics PUBLIC cxx_std_23)
target_link_libraries(hydrodynamics Eigen3::Eigen tinyxml2::tinyxml2)

set(EXAMPLES parse_model)

foreach(example IN ITEMS ${EXAMPLES})
    add_executable(${example} examples/${example}.cpp)
    add_dependencies(${example} hydrodynamics)
    target_link_libraries(${example} PUBLIC hydrodynamics)
    set_target_properties(
        ${example}
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/examples
    )
endforeach()

install(
    TARGETS hydrodynamics
    EXPORT export_hydrodynamics
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    FILE_SET HEADERS
)

ament_export_targets(export_hydrodynamics HAS_LIBRARY_TARGET)
ament_export_dependencies(Eigen3)

ament_package()
