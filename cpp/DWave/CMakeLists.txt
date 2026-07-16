cmake_minimum_required(VERSION 3.30)

project(DWave VERSION 1.0.0 LANGUAGES CXX)

list(APPEND CMAKE_PREFIX_PATH "$ENV{HOME}/usr/local")

find_package(OpenMP REQUIRED)
find_package(mrock 
    REQUIRED 
    COMPONENTS 
        symbolic_operators 
        iEoM 
        utility
)

add_executable(DWave 
    sources/DWave.cpp
)
target_sources(DWave PRIVATE
    sources/DWave/Model.cpp
    sources/DWave/T_C.cpp
)

# provided by the mrock library
mrock_generate_information_header(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/generated/dwave/info.h"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    OUT_HEADER MROCK_GENERATED_INFO_HEADER
    OUT_INCLUDE_DIR MROCK_GENERATED_INCLUDE_DIR
)
include(${CMAKE_CURRENT_LIST_DIR}/cmake/dwave_options.cmake)
# Sets include dir for the generated metadata header file
target_include_directories(dwave_options INTERFACE ${MROCK_GENERATED_INCLUDE_DIR})

target_link_libraries(DWave 
    PRIVATE  
        OpenMP::OpenMP_CXX 
        ${mrock_LIBRARIES}
        dwave_options
)

option(LATTICE_CUT_RESIDUALS "Compute residuals and eigenvectors of the transformed dynamical matrix." OFF)
if (LATTICE_CUT_RESIDUALS)
    target_compile_definitions(DWave PRIVATE LATTICE_CUT_RESIDUALS)
endif()
option(LATTICE_CUT_FULL_DIAG "Fully diagonalize the transformed dynamical matrix." OFF)
if (LATTICE_CUT_FULL_DIAG)
    target_compile_definitions(DWave PRIVATE LATTICE_CUT_FULL_DIAG)
endif()