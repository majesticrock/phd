# ModeCommutators/CMakeLists.txt

cmake_minimum_required(VERSION 3.30)

project(ModeCommutators VERSION 1.0)

#
# Set local settings and set the default build type as Release
#

list(APPEND CMAKE_PREFIX_PATH "$ENV{HOME}/usr/local")
set(Boost_DIR /usr/lib64/openmpi/lib/cmake/Boost-1.78.0)

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
message(STATUS "Build type is '${CMAKE_BUILD_TYPE}'.")

#
# Library specifics
#

add_executable(ModeCommutators 
    sources/ModeCommutators.cpp
)

target_compile_features(ModeCommutators
    INTERFACE
        cxx_std_20
)

find_package(mrock 
    REQUIRED 
    COMPONENTS 
        symbolic_operators 
)

target_link_libraries(ModeCommutators PUBLIC ${mrock_LIBRARIES})
