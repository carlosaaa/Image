cmake_minimum_required(VERSION 3.7)
project(Image)

set(CMAKE_CXX_STANDARD 98)


set(SOURCE_FILES src/EasyBMP/EasyBMP.cpp src/Log.cpp src/main.cpp)
INCLUDE_DIRECTORIES(src/EasyBMP)

add_executable(Image ${SOURCE_FILES})