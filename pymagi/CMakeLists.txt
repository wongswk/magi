cmake_minimum_required(VERSION 3.17)

set(pybind11_DIR ${CMAKE_CURRENT_LIST_DIR}/../include/pybind11/)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -flto -O3 -W -Wall -Wformat=2 -Wnull-dereference -Wold-style-cast -pedantic -fprofile-arcs -ftest-coverage -fPIC ")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER g++-6)

file(GLOB LIB_PYMAGI_HEADERS "src/*.h")
file(GLOB LIB_PYMAGI_SOURCE "src/*.cpp")

add_library(pymagi SHARED ${LIB_PYMAGI_SOURCE})
set_target_properties(pymagi PROPERTIES LINKER_LANGUAGE CXX)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(SUFFIX .cpython-36m-darwin.so)
    set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
else ()
    set(SUFFIX .cpython-36m-x86_64-linux-gnu.so)
endif ()

find_package (Python3 COMPONENTS NumPy)

target_include_directories(pymagi PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../cmagi/)
target_include_directories(pymagi SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include/)
target_include_directories(pymagi SYSTEM PUBLIC
        /usr/local/Cellar/python/3.6.4_4/Frameworks/Python.framework/Versions/3.6/include/python3.6m
        /usr/local/Cellar/python/3.6.5_1/Frameworks/Python.framework/Versions/3.6/include/python3.6m
        /usr/include/python3.6m/
        /usr/local/Frameworks/Python.framework/Headers
        )

target_include_directories(pymagi SYSTEM PUBLIC ${Python3_NumPy_INCLUDE_DIRS})


set_target_properties(pymagi PROPERTIES PREFIX "" SUFFIX "${SUFFIX}")
set_target_properties(pymagi PROPERTIES OUTPUT_NAME "pymagi")
find_library(CMAGI NAMES cmagi PATHS ${CMAKE_CURRENT_SOURCE_DIR}/../cmagi/)
target_link_libraries(pymagi ${CMAGI})

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../cmagi/)
