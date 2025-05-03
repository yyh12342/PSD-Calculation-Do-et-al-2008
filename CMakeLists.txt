cmake_minimum_required(VERSION 3.10)

project(PSDCalculator)

# C++17, no compiler extensions
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 실행 파일
add_executable(PSDCalculator
    src/main.cpp
    src/Potential.cpp
    src/SphereFinder.cpp
    src/Histogram.cpp)

# include 디렉터리
target_include_directories(PSDCalculator
    PRIVATE src)