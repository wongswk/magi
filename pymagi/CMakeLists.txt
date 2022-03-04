cmake_minimum_required(VERSION 3.16)

set(pybind11_DIR ${CMAKE_CURRENT_LIST_DIR}/../include/pybind11/)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -flto -O3 -W -Wall -Wformat=2 -Wnull-dereference -Wold-style-cast -pedantic -fprofile-arcs -ftest-coverage -fPIC ")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER g++)

file(GLOB LIB_PYMAGI_HEADERS "src/*.h")
file(GLOB LIB_PYMAGI_SOURCE "src/*.cpp")

add_library(pymagi SHARED ${LIB_PYMAGI_SOURCE})
set_target_properties(pymagi PROPERTIES LINKER_LANGUAGE CXX)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif ()

find_package (Python3 COMPONENTS NumPy Interpreter Development)

target_include_directories(pymagi PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../cmagi/)
target_include_directories(pymagi SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include/)
target_include_directories(pymagi SYSTEM PUBLIC ${Python3_INCLUDE_DIRS})
target_include_directories(pymagi SYSTEM PUBLIC ${Python3_NumPy_INCLUDE_DIRS})
set(SUFFIX ".${Python3_SOABI}.so")


set_target_properties(pymagi PROPERTIES PREFIX "" SUFFIX "${SUFFIX}")
set_target_properties(pymagi PROPERTIES OUTPUT_NAME "pymagi")
find_library(CMAGI NAMES cmagi PATHS ${CMAKE_CURRENT_SOURCE_DIR}/../cmagi/)
target_link_libraries(pymagi ${CMAGI})

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../cmagi/)
