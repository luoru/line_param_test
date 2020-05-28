cmake_minimum_required(VERSION 2.8)

project(line_param_test)

set(Eigen3_DIR "C:/Users/LUOZHONGFEI/ThirdParty/Eigen3/share/eigen3/cmake")

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(line_param_test line_param_test.cpp)