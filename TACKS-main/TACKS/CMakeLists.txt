cmake_minimum_required(VERSION 4.0)
project(proj__)

set(CMAKE_CXX_STANDARD 20)

include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(proj__ main.cpp
        include/Extractor.h
        include/divider.h
        include/benchmark.h
        include/BF.h
        include/k-means.h
        include/Parameter.h
        include/Tracedata.h
        SimpleLayer/Simple_Sketch.h
        SimpleLayer/benchmark3.h
        DSSketch/DSSketch.h
        DSSketch/benchmark4.h
        include/BF_2.h
        Log_layer/log_sketch.h
        "k-means++/km_sketch.h"
        k-means++/benchmark_km.h
)

target_compile_options(proj__ PRIVATE
        -O3
        -ffast-math
)