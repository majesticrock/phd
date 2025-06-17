cmake_minimum_required(VERSION 3.16)

project(test_boost_cmake VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# --- Dependencies ---
include(CMakeFindDependencyMacro)
find_dependency(Boost REQUIRED COMPONENTS serialization iostreams)

# --- Library Definition ---
add_library(test_boost_cmake INTERFACE)

target_include_directories(test_boost_cmake INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(test_boost_cmake INTERFACE 
    Boost::serialization
    Boost::iostreams
)
