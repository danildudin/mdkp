cmake_minimum_required(VERSION 3.17)
project(mdkp VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories("include/")

add_executable(mdkp src/main.cpp src/CMDArgs.cpp src/CORALSolver.cpp src/CoreData.cpp src/LPSolution.cpp src/Mdkp.cpp src/Problem.cpp src/utils.cpp)

target_link_libraries(mdkp glpk)
