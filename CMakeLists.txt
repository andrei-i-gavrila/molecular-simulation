cmake_minimum_required(VERSION 3.15)
project(lab1)

set(CMAKE_CXX_STANDARD 17)

set(NOT_OPTIMZED_SOURCES
        "1 MD Not Optimized/globaldata.c"
        "1 MD Not Optimized/globaldata.h"
        "1 MD Not Optimized/initializer.c"
        "1 MD Not Optimized/initializer.h"
        "1 MD Not Optimized/main.c"
        "1 MD Not Optimized/running.c"
        "1 MD Not Optimized/running.h"
        "1 MD Not Optimized/timing.c"
        "1 MD Not Optimized/timing.h"
        )

set(OPTIMIZED_SOURCES
        "1 MD Optimized with Verlet/globaldata.c"
        "1 MD Optimized with Verlet/globaldata.h"
        "1 MD Optimized with Verlet/initializer.c"
        "1 MD Optimized with Verlet/initializer.h"
        "1 MD Optimized with Verlet/main.c"
        "1 MD Optimized with Verlet/running.c"
        "1 MD Optimized with Verlet/running.h"
        "1 MD Optimized with Verlet/timing.c"
        "1 MD Optimized with Verlet/timing.h")


add_executable(non_optimized ${NOT_OPTIMZED_SOURCES})
add_executable(optimized ${OPTIMIZED_SOURCES})

