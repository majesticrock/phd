# Lindhard/CMakeLists.txt

cmake_minimum_required(VERSION 3.30)

project(Lindhard VERSION 1.0)

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

find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED COMPONENTS serialization)
find_package(mrock 
    REQUIRED 
    COMPONENTS 
        utility
)

add_executable(Lindhard 
    sources/Lindhard.cpp
)

target_compile_features(Lindhard
    INTERFACE
        cxx_std_20
)

target_link_libraries(Lindhard PUBLIC OpenMP::OpenMP_CXX ${mrock_LIBRARIES})
