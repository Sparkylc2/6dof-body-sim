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

# --- clarabel (conic solver, used by scvx) ---
# the solver is written in rust, Clarabel.cpp is the c++/eigen wrapper around it.
# its own cmake reinstalls cbindgen and builds the examples every build, so this only fetches the
# source and runs cargo on the wrapper crate directly (the headers it ships are already generated)
include(FetchContent)
FetchContent_Declare(clarabel
  GIT_REPOSITORY https://github.com/oxfordcontrol/Clarabel.cpp.git
  GIT_TAG v0.11.1
  GIT_SHALLOW TRUE
  SOURCE_SUBDIR no-cmake # no CMakeLists here, so it's fetched but not add_subdirectory'd
)
FetchContent_MakeAvailable(clarabel)

find_program(CARGO cargo REQUIRED HINTS $ENV{HOME}/.cargo/bin)
set(CLARABEL_CARGO_DIR ${CMAKE_BINARY_DIR}/clarabel_cargo)
# always a release build, the rust code isn't something being debugged here and debug is very slow
set(CLARABEL_LIB ${CLARABEL_CARGO_DIR}/release/libclarabel_c.a)

add_custom_command(
  OUTPUT ${CLARABEL_LIB}
  COMMAND ${CARGO} build --release --target-dir ${CLARABEL_CARGO_DIR}
  WORKING_DIRECTORY ${clarabel_SOURCE_DIR}/rust_wrapper
  COMMENT "building clarabel (rust, first build takes a minute)"
)
add_custom_target(clarabel_cargo DEPENDS ${CLARABEL_LIB})

add_library(clarabel STATIC IMPORTED GLOBAL)
set_target_properties(clarabel PROPERTIES IMPORTED_LOCATION ${CLARABEL_LIB})
# imported include dirs are SYSTEM, so its headers don't trip -Werror
target_include_directories(clarabel INTERFACE ${clarabel_SOURCE_DIR}/include)
target_link_libraries(clarabel INTERFACE Eigen3::Eigen)
add_dependencies(clarabel clarabel_cargo)


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
  src/control/blocks.cpp
  src/control/autopilots/autopilot.cpp
  src/control/autopilots/lateral_autopilot.cpp
  src/control/autopilots/attitude_control.cpp
  src/control/autopilots/longitudinal_autopilot.cpp
  src/world.cpp
  src/control/solvers/conic.cpp
  src/control/solvers/problem_builder.cpp
)
target_include_directories(bsim_core PUBLIC include)
target_link_libraries(bsim_core PUBLIC Eigen3::Eigen PRIVATE raylib clarabel)
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
  src/scenes/scvx_scene.cpp
)
target_link_libraries(bsim PRIVATE bsim_render)
target_compile_options(bsim PRIVATE ${BSIM_WARNINGS})

enable_testing()
add_executable(bsim_tests tests/test_main.cpp)
target_link_libraries(bsim_tests PRIVATE bsim_core)
target_compile_options(bsim_tests PRIVATE ${BSIM_WARNINGS})
add_test(NAME bsim_tests COMMAND bsim_tests)
