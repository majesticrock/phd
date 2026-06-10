# Lindhard/CMakeLists.txt

cmake_minimum_required(VERSION 3.16)

project(Lindhard VERSION 1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_PREFIX_PATH "~/usr/local")

include_directories(
    $ENV{HOME}/usr/local/include
)

find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED COMPONENTS serialization)
find_package(mrock REQUIRED)

add_executable(Lindhard 
    sources/Lindhard.cpp
)

# Check if using GCC and adjust C++ standard flag if necessary
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    # Check GCC version
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11.0)
        message(FATAL_ERROR "GCC version ${CMAKE_CXX_COMPILER_VERSION} is not supported. GCC 12.1 or newer is required for C++20 support.")
    endif()
    # Use the appropriate flag for GCC
    target_compile_options(Lindhard PRIVATE
        -Wall
        -Wno-sign-compare
        -fopenmp
        -march=native
        -O3
    )
else()
    message(FATAL_ERROR "Unsupported compiler ${CMAKE_CXX_COMPILER_ID}. Only GCC is supported.")
endif()


target_include_directories(Lindhard PRIVATE ${Boost_INCLUDE_DIRS} ${mrock_INCLUDE_DIRS})
target_link_libraries(Lindhard PUBLIC OpenMP::OpenMP_CXX ${mrock_LIBRARIES})
