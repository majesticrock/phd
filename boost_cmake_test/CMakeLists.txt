cmake_minimum_required(VERSION 3.16)

project(test_boost_cmake VERSION 1.0.0 LANGUAGES CXX)

set(Boost_DIR /usr/lib64/openmpi/lib/cmake/Boost-1.78.0)# the cluster needs help finding the correct version of boost

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# --- Dependencies ---
include(CMakeFindDependencyMacro)
find_dependency(Boost REQUIRED COMPONENTS serialization iostreams)

# --- Library Definition ---
add_executable(test_boost_cmake test.cpp)

target_link_libraries(test_boost_cmake PRIVATE 
    Boost::serialization
    Boost::iostreams
)
target_include_directories(test_boost_cmake PRIVATE ${Boost_INCLUDE_DIRS})

message(WARNING "Loading boost libraries from ${Boost_INCLUDE_DIRS}")