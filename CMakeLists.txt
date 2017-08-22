cmake_minimum_required(VERSION 3.8)
project(Neelix)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp AuxiliaryField.cpp AuxiliaryField.h)
add_executable(Neelix ${SOURCE_FILES} AuxiliaryField.cpp AuxiliaryField.h)