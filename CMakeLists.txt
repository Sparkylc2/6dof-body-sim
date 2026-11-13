cmake_minimum_required(VERSION 3.24)
project(bsim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

add_compile_options($<$<CONFIG:Debug>:-fsanitize=address,undefined>)
add_link_options($<$<CONFIG:Debug>:-fsanitize=address,undefined>)
add_compile_definitions($<$<CONFIG:Debug>:EIGEN_INITIALIZE_MATRICES_BY_NAN>)

set(BSIM_WARNINGS -Wall -Wextra -Wpedantic -Werror)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(raylib REQUIRED)


add_library(bsim_core STATIC
  src/frames.cpp
  src/force_models.cpp
  src/force_models/fixed_wing_model.cpp
  src/force_models/spacecraft_model.cpp
  src/body.cpp
  src/trim.cpp
  src/spacecraft.cpp
  src/body_system.cpp
  src/air_system.cpp
  src/sensors.cpp
  src/control.cpp
  src/autopilots/autopilot.cpp
  src/autopilots/lateral_autopilot.cpp
  src/autopilots/attitude_control.cpp
  src/world.cpp
)
target_include_directories(bsim_core PUBLIC include)
target_link_libraries(bsim_core PUBLIC Eigen3::Eigen PRIVATE raylib)
target_compile_options(bsim_core PRIVATE ${BSIM_WARNINGS})

add_library(bsim_render STATIC
  src/render/render.cpp
)
target_include_directories(bsim_render PUBLIC src)
target_link_libraries(bsim_render PUBLIC bsim_core raylib)
target_compile_options(bsim_render PRIVATE ${BSIM_WARNINGS})

add_executable(bsim
  src/main.cpp
  src/scenes/plane_scene.cpp
  src/scenes/spacecraft_scene.cpp
  src/scenes/sindy_scene.cpp
)
target_link_libraries(bsim PRIVATE bsim_render)
target_compile_options(bsim PRIVATE ${BSIM_WARNINGS})

enable_testing()
add_executable(bsim_tests tests/test_main.cpp)
target_link_libraries(bsim_tests PRIVATE bsim_core)
target_compile_options(bsim_tests PRIVATE ${BSIM_WARNINGS})
add_test(NAME bsim_tests COMMAND bsim_tests)
