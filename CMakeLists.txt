cmake_minimum_required(VERSION 3.9)
project(gpds_multi_lang)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_COMPILER g++-6)
set(CMAKE_MAKE_PROGRAM)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/gpds_cpp)
#add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/pygpds)
