cmake_minimum_required (VERSION 3.0)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fext-numeric-literals")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${EXT_PROJECTS_DIR}/bin)

include_directories(${MCMC_SOURCE_DIR}/include)

## Compiler flags
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")        ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
endif()

add_library(Dynamics Dynamics.cpp)

add_library(FourierSeries FourierSeries.cpp)

add_library(Grid Grid.cpp)
target_link_libraries( Grid FourierSeries )

add_library(LangevinDynamics LangevinDynamics.cpp)
target_link_libraries( LangevinDynamics FourierSeries Dynamics Options)

add_library(CurvedSurfaceDynamics CurvedSurfaceDynamics.cpp)
target_link_libraries( CurvedSurfaceDynamics FourierSeries Dynamics Options)

add_library(BridgeDynamics BridgeDynamics.cpp)
target_link_libraries( BridgeDynamics FourierSeries Dynamics Options )

add_library(LikelihoodFreeMCMC LikelihoodFreeMCMC.cpp)

add_library(Options Options.cpp)

add_library(OUDynamics OUDynamics.cpp)
target_link_libraries(OUDynamics Dynamics)