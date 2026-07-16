# Lindhard/CMakeLists.txt

cmake_minimum_required(VERSION 3.30)

project(Lindhard VERSION 1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_PREFIX_PATH "$ENV{HOME}/usr/local")

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
