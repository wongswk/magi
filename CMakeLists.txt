cmake_minimum_required(VERSION 3.9)
project(magi_multi_lang)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_COMPILER g++-6)
set(CMAKE_MAKE_PROGRAM)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(/Users/shihaoyang/R/library/Rcpp/include)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/cmagi)
#add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/pymagi)
