# HubbardMeanField/CMakeLists.txt

cmake_minimum_required(VERSION 3.19)

project(HubbardMeanField VERSION 1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_PREFIX_PATH "~/usr/local")

if(CLUSTER_BUILD)
    include(${CMAKE_SOURCE_DIR}/cmake/ClusterCompilerFlags.cmake)
else()
    include(${CMAKE_SOURCE_DIR}/cmake/DefaultCompilerFlags.cmake)
endif()

include_directories(
    $ENV{HOME}/usr/local/include
)

add_executable(HubbardMeanField 
    sources/HubbardMeanField.cpp
)

SET_COMPILER_FLAGS(HubbardMeanField)

find_package(MPI REQUIRED)
target_link_libraries(HubbardMeanField PRIVATE MPI::MPI_CXX)

find_package(Boost REQUIRED COMPONENTS serialization iostreams)
target_include_directories(HubbardMeanField PRIVATE ${Boost_INCLUDE_DIRS})

# System libraries
find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Utility REQUIRED)

# SymbolicOperators library
find_package(SymbolicOperators REQUIRED CONFIG)
target_include_directories(HubbardMeanField PRIVATE ${SymbolicOperators_INCLUDE_DIRS})
target_link_libraries(HubbardMeanField PRIVATE SymbolicOperators::SymbolicOperators)

# Utility library
find_package(Utility REQUIRED)
target_include_directories(HubbardMeanField PRIVATE ${Utility_INCLUDE_DIRS})
target_link_libraries(HubbardMeanField PRIVATE Utility::Utility)

# subdirs
add_subdirectory(sources/Hubbard/Helper)
add_subdirectory(sources/Handler)
add_subdirectory(sources/Hubbard/Models)
add_subdirectory(sources/Hubbard/DensityOfStates)

target_link_libraries(HubbardMeanField PRIVATE
    ${Boost_LIBRARIES}
    OpenMP::OpenMP_CXX
    ZLIB::ZLIB
    Models
    DensityOfStates
    Handler
    Helper
)

target_link_options(HubbardMeanField PRIVATE "-Wl,--trace-symbol")