# the matlab bridge (only built with -DBSIM_MATLAB=ON, from the root CMakeLists)
#
#   bsim_mex     the mex the +bsim package calls (self contained)
#   libbsim_c    the c api as a shared library, for simulink's C Caller / C Function blocks
#
# both land in matlab/bin so setup.m only has to add that one folder

# the debug build has address sanitizer on everything, and loading that into matlab crashes it
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  message(FATAL_ERROR
    "the matlab bridge can't be a Debug build (asan crashes matlab), configure it with "
    "-DCMAKE_BUILD_TYPE=Release (e.g. in build-release)")
endif()

# finds the matlab on the PATH, or set Matlab_ROOT_DIR to pick one
find_package(Matlab REQUIRED COMPONENTS MX_LIBRARY)

set(BSIM_MATLAB_BIN ${PROJECT_SOURCE_DIR}/matlab/bin)

# the shared code, compiled once and put into both
add_library(bsim_bridge OBJECT bsim_c.cpp)
target_include_directories(bsim_bridge PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(bsim_bridge PUBLIC bsim_models)
target_compile_options(bsim_bridge PRIVATE ${BSIM_WARNINGS})
set_target_properties(bsim_bridge PROPERTIES POSITION_INDEPENDENT_CODE ON)

# --- mex ---
# R2018a is the interleaved complex api (mxGetDoubles)
# (only the c api libraries, it doesn't need the engine or data array ones)
matlab_add_mex(NAME bsim_mex SRC bsim_mex.cpp R2018a NO_IMPLICIT_LINK_TO_MATLAB_LIBRARIES
  LINK_TO bsim_bridge ${Matlab_MEX_LIBRARY} ${Matlab_MX_LIBRARY})
# matlab's headers don't get along with -Wpedantic, so the mex only gets the main warnings
target_compile_options(bsim_mex PRIVATE -Wall -Wextra -Werror)
set_target_properties(bsim_mex PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${BSIM_MATLAB_BIN}
  RUNTIME_OUTPUT_DIRECTORY ${BSIM_MATLAB_BIN})
if(APPLE)
  # matlab_add_mex always adds the c++ mex api's export list too, which the newer apple linker
  # refuses for a c api mex (those symbols don't exist), so this only exports what a c api mex has
  set_target_properties(bsim_mex PROPERTIES LINK_FLAGS
    "-Wl,-exported_symbols_list,${Matlab_EXTERN_LIBRARY_DIR}/c_exportsmexfileversion.map")
endif()

# --- shared library for simulink ---
# simulink links its own simulation target against this, so the install name is the absolute path
# of the built file, which lets it load from wherever simulink builds (rebuild if the repo moves)
add_library(bsim_c SHARED)
target_link_libraries(bsim_c PRIVATE bsim_bridge)
set_target_properties(bsim_c PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${BSIM_MATLAB_BIN}
  BUILD_WITH_INSTALL_NAME_DIR ON
  INSTALL_NAME_DIR ${BSIM_MATLAB_BIN})
