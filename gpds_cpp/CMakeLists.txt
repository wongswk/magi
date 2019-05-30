cmake_minimum_required(VERSION 3.9)
project(gpds_cpp)

set(CMAKE_CXX_STANDARD 11)


if (IS_DIRECTORY $ENV{HOME}/R_LIBS)
    set(R_LIBS_USER $ENV{HOME}/R_LIBS)
else ()
    set(R_LIBS_USER $ENV{HOME}/Library/R/3.4/library/)
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
find_library(BLAS NAMES openblas /usr/local/opt/openblas/lib)
find_library(LAPACK NAMES lapack PATHS /usr/lib/)
find_library(RBLAS NAMES Rblas PATHS ${R_PATH}/Resources/lib/ ${R_PATH}/lib/)
find_library(RLAPACK NAMES Rlapack PATHS ${R_PATH}/Resources/lib/ ${R_PATH}/lib/)

include_directories(/usr/local/opt/openblas/include)

add_library(gpds_cpp SHARED library.cpp library.h hmc.cpp hmc.h classDefinition.cpp classDefinition.h)  #  band.cpp band.h
target_include_directories(gpds_cpp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../include /usr/local/opt/openblas/include)

SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath=../lib,-rpath=/usr/local/opt/openblas/lib")

