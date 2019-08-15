
set(pybind11_DIR ${CMAKE_CURRENT_LIST_DIR}/../include/pybind11/)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -flto -O3 -W -Wall -Werror -Wformat=2 -Wnull-dereference -Wold-style-cast -pedantic -fprofile-arcs -ftest-coverage -fPIC ")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_COMPILER g++-6)

file(GLOB LIB_PYGPDS_HEADERS "src/*.h")
file(GLOB LIB_PYGPDS_SOURCE "src/*.cpp")

add_library(pygpds SHARED ${LIB_PYGPDS_SOURCE})
set_target_properties(pygpds PROPERTIES LINKER_LANGUAGE CXX)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(SUFFIX .cpython-36m-darwin.so)
    set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
else ()
    set(SUFFIX .cpython-36m-x86_64-linux-gnu.so)
endif ()

target_include_directories(pygpds PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../gpds_cpp/)
target_include_directories(pygpds SYSTEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include/)
target_include_directories(pygpds SYSTEM PUBLIC
        /usr/local/Cellar/python/3.6.4_4/Frameworks/Python.framework/Versions/3.6/include/python3.6m
        /usr/local/Cellar/python/3.6.5_1/Frameworks/Python.framework/Versions/3.6/include/python3.6m
        /usr/include/python3.6m/)

add_dependencies(pygpds cgpds)

set_target_properties(pygpds PROPERTIES PREFIX "" SUFFIX "${SUFFIX}")
set_target_properties(pygpds PROPERTIES OUTPUT_NAME "pygpds")
find_library(CGPDS NAMES cgpds PATHS ${CMAKE_CURRENT_SOURCE_DIR}/../gpds_cpp/)
target_link_libraries(pygpds ${CGPDS})

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../gpds_cpp/)
