cmake_minimum_required(VERSION 3.10)

project(PoreSizeDistributionCalculator)

set(CMAKE_CXX_STANDARD 17)

add_executable(PoreSizeDistributionCalculator 
    src/main.cpp 
    Geometry.h # 우선순위 때문에 추가
    # src/PDBParser.h
    # src/PDBParser.cpp 
    src/Potential.h
    src/Potential.cpp 
    src/SphereFinder.h
    src/SphereFinder.cpp 
    src/Histogram.h
    src/Histogram.cpp)