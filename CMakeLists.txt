cmake_minimum_required(VERSION 3.8)
project(Neelix)

set(CMAKE_CXX_STANDARD 11)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

set(SOURCE_FILES main.cpp AuxiliaryField.cpp AuxiliaryField.h FermionMatrix.cpp FermionMatrix.h)
add_executable(Neelix ${SOURCE_FILES} AuxiliaryField.cpp AuxiliaryField.h FermionMatrix.cpp FermionMatrix.h)
target_link_libraries(Neelix ${ARMADILLO_LIBRARIES})