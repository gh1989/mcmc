set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fext-numeric-literals")
include_directories( ${MCMC_SOURCE_DIR}/include )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${EXT_PROJECTS_DIR}/bin)

add_executable(LFOrnsteinUhlenbeck LFOrnsteinUhlenbeck.cpp)
target_link_libraries( LFOrnsteinUhlenbeck gsl gslcblas Algorithm LikelihoodFreeScheme Options OUPathScheme OUParams FourierSeries)

add_executable(LFLangevin LFLangevin.cpp)
target_link_libraries( LFLangevin gsl gslcblas Algorithm LikelihoodFreeScheme Options LangevinPathScheme LangevinParams FourierSeries)
