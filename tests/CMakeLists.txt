set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fext-numeric-literals")
include_directories( ${MCMC_SOURCE_DIR}/include )

add_library(fftw3 STATIC IMPORTED)
set_property(TARGET fftw3 PROPERTY IMPORTED_LOCATION /usr/local/lib/libfftw3.a)
add_library(gsl STATIC IMPORTED)
set_property(TARGET gsl PROPERTY IMPORTED_LOCATION /usr/local/lib/libgsl.a)
add_library(gslcblas STATIC IMPORTED)
set_property(TARGET gslcblas PROPERTY IMPORTED_LOCATION /usr/local/lib/libgslcblas.a)

add_executable(TestNumericLimits TestNumericLimits.cpp)

add_executable(TestFourierSeries TestFourierSeries.cpp)
target_link_libraries(TestFourierSeries FourierSeries)

add_executable(TestOptions TestOptions.cpp)
target_link_libraries(TestOptions Options)

add_executable(TestLangevin TestLangevin.cpp)
target_link_libraries(TestLangevin LangevinDynamics gsl gslcblas Options)

add_executable(TestRandom TestRandom.cpp)
target_link_libraries(TestRandom gsl gslcblas)

add_executable(TestVector TestVector.cpp)
target_link_libraries(TestVector LangevinDynamics BridgeDynamics gsl gslcblas Options)

add_executable(TestDiffusionBridge TestDiffusionBridge.cpp)
target_link_libraries(TestDiffusionBridge LangevinDynamics BridgeDynamics gsl gslcblas Options)