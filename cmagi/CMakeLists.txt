cmake_minimum_required(VERSION 3.9)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_COMPILER g++-6)


if (IS_DIRECTORY $ENV{HOME}/R/library)
    set(R_LIBS_USER $ENV{HOME}/R/library)
else ()
    set(R_LIBS_USER $ENV{HOME}/Library/R/3.6/library/)
endif ()


if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(R_PATH /Library/Frameworks/R.framework)
else ()

    # add fopenmp support only to linux
    # NOTE: shadow warnings disabled on linux since gcc will throw errors about constructor syntax
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -fopenmp -Wno-shadow")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fopenmp -Wno-shadow")

    # AWS
    if (IS_DIRECTORY /usr/lib64/R/)
        set(R_PATH /usr/lib64/R/)
    elseif (IS_DIRECTORY /usr/lib/R)
        set(R_PATH /usr/lib/R)
    endif ()
endif ()


find_library(ARMADILLO NAMES armadillo PATHS ${CMAKE_CURRENT_SOURCE_DIR}/../lib/)
find_library(BLAS NAMES blas PATHS /usr/lib/)
find_library(LAPACK NAMES lapack PATHS /usr/lib/)
find_library(RBLAS NAMES Rblas PATHS ${R_PATH}/Resources/lib/ ${R_PATH}/lib/)
find_library(RLAPACK NAMES Rlapack PATHS ${R_PATH}/Resources/lib/ ${R_PATH}/lib/)
find_library(PTHREAD NAMES pthread PATHS /usr/lib/)


file(GLOB LIB_MAGI_HEADERS "*.h")
file(GLOB LIB_MAGI_SOURCE "*.cpp")

add_library(cmagi SHARED ${LIB_MAGI_HEADERS} ${LIB_MAGI_SOURCE})
target_include_directories(cmagi PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
#add_definitions(-DARMA_BLAS_LONG_LONG)
set_target_properties(cmagi PROPERTIES COMPILE_DEFINITIONS "ARMA_BLAS_LONG_LONG")
target_link_libraries(cmagi ${BLAS} ${ARMADILLO} ${LAPACK} ${PTHREAD})

SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
add_executable(magi_main ${CMAKE_CURRENT_SOURCE_DIR}/tests/main.cpp)
target_include_directories(magi_main PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include)
target_link_libraries(magi_main cmagi ${ARMADILLO})
