# ContinuumSystem/CMakeLists.txt

cmake_minimum_required(VERSION 3.19)

project(ContinuumSystem VERSION 1.0)

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

add_executable(ContinuumSystem 
    sources/ContinuumSystem.cpp
)
set(ADD_SOURCES
    sources/Continuum/ModeHelper.cpp
    sources/Continuum/SCModel.cpp
    sources/Continuum/MomentumRanges.cpp
    sources/Continuum/SplineContainer.cpp
    sources/Continuum/ModelInitializer.cpp
    sources/Continuum/Incrementer.cpp
)
target_sources(ContinuumSystem PRIVATE ${ADD_SOURCES})

SET_COMPILER_FLAGS(ContinuumSystem)

find_package(MPI REQUIRED)
target_link_libraries(ContinuumSystem PRIVATE MPI::MPI_CXX)

find_package(Boost REQUIRED COMPONENTS serialization iostreams)
target_include_directories(ContinuumSystem PRIVATE ${Boost_INCLUDE_DIRS})

find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)

find_package(SymbolicOperators REQUIRED CONFIG)
target_include_directories(ContinuumSystem PRIVATE ${SymbolicOperators_INCLUDE_DIRS})
target_link_libraries(ContinuumSystem PRIVATE SymbolicOperators::SymbolicOperators)

find_package(Utility REQUIRED)
target_include_directories(ContinuumSystem PRIVATE ${Utility_INCLUDE_DIRS})
target_link_libraries(ContinuumSystem PRIVATE Utility::Utility)

target_link_libraries(ContinuumSystem PRIVATE OpenMP::OpenMP_CXX ${Boost_LIBRARIES} ZLIB::ZLIB)

