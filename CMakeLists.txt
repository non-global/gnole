cmake_minimum_required(VERSION 2.8.9)
project(ngl)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin) 

include_directories(${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src)

set(SOURCES
	${PROJECT_SOURCE_DIR}/src/CmdLine.cc
	${PROJECT_SOURCE_DIR}/src/Event.cc
	${PROJECT_SOURCE_DIR}/src/Momentum.cc
	${PROJECT_SOURCE_DIR}/src/Parameters.cc
	${PROJECT_SOURCE_DIR}/src/Shower.cc
	${PROJECT_SOURCE_DIR}/src/helpers.cc)
file(GLOB HEADERS ${PROJECT_SOURCE_DIR}/include/*.hh)
add_executable(gnole
        ${PROJECT_SOURCE_DIR}/src/gnole.cc
	${SOURCES} ${HEADERS})
add_executable(test_evolgrid
        ${PROJECT_SOURCE_DIR}/src/test_evolgrid.cc
	${SOURCES} ${HEADERS})

# just for example add some compiler flags
add_definitions(-Wall -g -O2)

find_package(GSL REQUIRED)
option(NNET "Use neural network implementation" OFF)
if(NNET)
find_package(Torch REQUIRED)
target_link_libraries(gnole ${TORCH_LIBRARIES} GSL::gsl GSL::gslcblas)
target_link_libraries(test_evolgrid ${TORCH_LIBRARIES} GSL::gsl GSL::gslcblas)
add_definitions(-DNNET)
else()
target_link_libraries(gnole GSL::gsl GSL::gslcblas)
target_link_libraries(test_evolgrid GSL::gsl GSL::gslcblas)
endif(NNET)
unset(NNET)
 
