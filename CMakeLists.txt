cmake_minimum_required(VERSION 3.18)

project(atom_search_cpp VERSION 1.0 LANGUAGES CXX)

# 包含检测编译器支持的模块
include(CheckCXXCompilerFlag)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")

# gprof
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
#set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
# 设置 Python 的头文件路径
#set(Python3_INCLUDE_DIRS "/home/public/miniconda3/envs/develop/include/python3.8")
# 设置 Python 的库路径
#set(Python3_LIBRARIES "/home/public/miniconda3/envs/develop/lib/libpython3.8.so")

find_package(Python3 REQUIRED COMPONENTS Development)
find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)
add_definitions(-DPYTHON_SCRIPT_PATH="${CMAKE_SOURCE_DIR}/script")

if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# 设置 EIGEN_USE_OPENMP 宏
#add_definitions(-DEIGEN_USE_OPENMP)
message(STATUS "Python3 Include Dirs: ${Python3_INCLUDE_DIRS}")
message(STATUS "Python3 Libraries: ${Python3_LIBRARIES}")

# Python 头文件
#include_directories("/home/public/miniconda3/envs/develop/include/python3.8")
# Eigen3
#include_directories("/usr/local/include/eigen3")
# Python 库文件
#link_directories("/home/public/miniconda3/envs/develop/lib/python3.8/site-packages")
#link_libraries("/home/public/miniconda3/envs/develop/lib/libpython3.8.so")

add_executable(atom-search-cpp main.cpp atom.cpp chemio.cpp chemtable.cpp curve.cpp graph.cpp grid.cpp mathfunc.cpp spreading.cpp optimizer.cpp loadlib.cpp crosslinker.cpp)
target_include_directories(atom-search-cpp PRIVATE ${Python3_INCLUDE_DIRS} Eigen3::Eigen)
target_link_libraries(atom-search-cpp PRIVATE ${Python3_LIBRARIES} OpenMP::OpenMP_CXX dl Eigen3::Eigen)

# 现在设置这个目标的输出名称
set_target_properties(atom-search-cpp PROPERTIES OUTPUT_NAME "polychef")