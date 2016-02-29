cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(ProcData)

######## locate ROOT install, headers, includes
include(FindROOT.cmake)

set(PG_LIB $ENV{PG4_LIB_DIR})
set(PG_INC "$ENV{PG4_CODE}/include/Output")

include_directories(${ROOT_INCLUDE_DIR} ${PG_INC})
link_directories(${ROOT_LIBRARY_DIR} ${PG_LIB})
add_compile_options(-std=c++11)
add_executable(ProcData ProcData.cpp)
target_link_libraries(ProcData -lEventLib)
target_link_libraries(ProcData ${ROOT_LIBRARIES})
MESSAGE( STATUS "DIC " ${ROOT_INCLUDE_DIR})
#set(CMAKE_BUILD_TYPE release)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_BUILD_TYPE debug)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3 -Wextra -Wpedantic")
