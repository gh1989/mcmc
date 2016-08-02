set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fext-numeric-literals")
include_directories( ${MCMC_SOURCE_DIR}/include )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${EXT_PROJECTS_DIR}/bin)

#add_library(fftw3 STATIC IMPORTED /usr/local/lib/libfftw3.a)
#add_library(gsl STATIC IMPORTED /usr/local/lib/libgsl.a)
#add_library(gslcblas STATIC IMPORTED /usr/local/lib/libgslcblas.a)

add_library(fftw3 STATIC IMPORTED)
set_property(TARGET fftw3 PROPERTY IMPORTED_LOCATION /usr/local/lib/libfftw3.a)
add_library(gsl STATIC IMPORTED)
set_property(TARGET gsl PROPERTY IMPORTED_LOCATION /usr/local/lib/libgsl.a)
add_library(gslcblas STATIC IMPORTED)
set_property(TARGET gslcblas PROPERTY IMPORTED_LOCATION /usr/local/lib/libgslcblas.a)

add_executable(LFOrnsteinUhlenbeck LFOrnsteinUhlenbeck.cpp)
target_link_libraries( LFOrnsteinUhlenbeck gsl gslcblas LikelihoodFreeMCMC Options OUDynamics FourierSeries)

add_executable(LFLangevin LFLangevin.cpp)
target_link_libraries( LFLangevin gsl gslcblas LikelihoodFreeMCMC Options LangevinDynamics FourierSeries)

add_executable(SequentialLFMCMCLangevin SequentialLFMCMC.cpp)
target_link_libraries(SequentialLFMCMCLangevin gsl gslcblas SequentialLFMCMC Options LangevinDynamics FourierSeries)
