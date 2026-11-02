cmake_minimum_required(VERSION 3.24)
project(uav_sim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

add_compile_options($<$<CONFIG:Debug>:-fsanitize=address,undefined>)
add_link_options($<$<CONFIG:Debug>:-fsanitize=address,undefined>)
add_compile_definitions($<$<CONFIG:Debug>:EIGEN_INITIALIZE_MATRICES_BY_NAN>)

set(UAV_WARNINGS -Wall -Wextra -Wpedantic -Werror)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(raylib REQUIRED)


add_library(uav_core STATIC
  src/frames.cpp
  src/force_models.cpp
  src/force_models/fixed_wing_model.cpp
  src/body.cpp
  src/body_system.cpp
  src/air_system.cpp
  src/sensors.cpp
  src/control.cpp
  src/autopilot.cpp
  src/world.cpp
)
target_include_directories(uav_core PUBLIC include)
target_link_libraries(uav_core PUBLIC Eigen3::Eigen PRIVATE raylib)
target_compile_options(uav_core PRIVATE ${UAV_WARNINGS})

add_library(uav_render STATIC
  src/render/render.cpp
)
target_include_directories(uav_render PUBLIC src)
target_link_libraries(uav_render PUBLIC uav_core raylib)
target_compile_options(uav_render PRIVATE ${UAV_WARNINGS})

add_executable(uav_sim src/main.cpp)
target_link_libraries(uav_sim PRIVATE uav_render)
target_compile_options(uav_sim PRIVATE ${UAV_WARNINGS})

enable_testing()
add_executable(uav_tests tests/test_main.cpp)
target_link_libraries(uav_tests PRIVATE uav_core)
target_compile_options(uav_tests PRIVATE ${UAV_WARNINGS})
add_test(NAME uav_tests COMMAND uav_tests)
