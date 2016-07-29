cmake_minimum_required (VERSION 3.0)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fext-numeric-literals")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${EXT_PROJECTS_DIR}/bin)

include_directories(${MCMC_SOURCE_DIR}/include)

add_library(fftw3 STATIC IMPORTED /usr/local/lib/libfftw3.a)
add_library(gsl STATIC IMPORTED /usr/local/lib/libgsl.a)
add_library(gslcblas STATIC IMPORTED /usr/local/lib/libgslcblas.a)

add_library(Algorithm Algorithm.cpp)

add_library(FourierSeries FourierSeries.cpp)

add_library(Grid Grid.cpp)
target_link_libraries( Grid FourierSeries )

add_library(LangevinParams LangevinParams.cpp)
target_link_libraries( LangevinParams FourierSeries gsl gslcblas )

add_library(LangevinPathScheme LangevinPathScheme.cpp)

add_library(LikelihoodFreeScheme LikelihoodFreeScheme.cpp)

add_library(Options Options.cpp)

add_library(OUParams OUParams.cpp)

add_library(OUPathScheme OUPathScheme.cpp)

