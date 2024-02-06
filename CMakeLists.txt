cmake_minimum_required(VERSION 3.0.0)
project(buildin VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()

add_executable(buildin lab2_1.cpp test.cpp Matrix.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
