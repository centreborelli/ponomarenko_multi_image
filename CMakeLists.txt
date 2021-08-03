cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_STANDARD 14)
project(ponomarenko_multi_image)


# debug or release
SET(CMAKE_BUILD_TYPE "Release")
# SET(CMAKE_BUILD_TYPE "Debug")

SET(PLATFORM_COMPILE_FLAGS "-msse -msse2 -msse3 -msse4 -mssse3")  # SSE only

SET(COMMON_COMPILE_FLAGS "-fopenmp -Wall -Wextra -ansi")

SET(CMAKE_CXX_FLAGS "${PLATFORM_COMPILE_FLAGS} ${COMMON_COMPILE_FLAGS} ${MVS_FLAGS} $ENV{CXXFLAGS}")
SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -g2 -ggdb")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3")

# include MIPP headers
include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}/framework
    ${CMAKE_CURRENT_SOURCE_DIR}
)

add_library(framework
    framework/CFramework.cpp
    framework/CHistogram.cpp
    framework/CImage.cpp
    framework/CInspector.cpp
    framework/framework.cpp
    framework/libparser.cpp
    framework/operations.cpp
)

add_library(curve_filter
    curve_filter.cpp
)

add_executable(ponomarenko_multi
    algo_multi.cpp
)
target_link_libraries(ponomarenko_multi
    framework curve_filter
    png 
    m 
    fftw3f 
    fftw3f_threads 
    gomp
)

add_executable(ponomarenko_single
    algo_single.cpp
)
target_link_libraries(ponomarenko_single
    framework curve_filter
    png 
    m 
    fftw3f 
    fftw3f_threads 
    gomp
)

