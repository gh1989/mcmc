cmake_minimum_required (VERSION 3.0)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fext-numeric-literals")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${EXT_PROJECTS_DIR}/bin)

add_library(fftw3 STATIC IMPORTED)
set_property(TARGET fftw3 PROPERTY IMPORTED_LOCATION /usr/local/lib/libfftw3.a)
add_library(gsl STATIC IMPORTED)
set_property(TARGET gsl PROPERTY IMPORTED_LOCATION /usr/local/lib/libgsl.a)
add_library(gslcblas STATIC IMPORTED)
set_property(TARGET gslcblas PROPERTY IMPORTED_LOCATION /usr/local/lib/libgslcblas.a)

include_directories(${MCMC_SOURCE_DIR}/include)

add_executable(LFMCMCTuning LFMCMCTuning.cpp)
target_link_libraries( LFMCMCTuning gsl gslcblas LikelihoodFreeMCMC Options LangevinDynamics FourierSeries)