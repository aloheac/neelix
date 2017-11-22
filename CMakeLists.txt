cmake_minimum_required(VERSION 3.8)
project(Neelix)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

set(SOURCE_FILES main.cpp AuxiliaryField.cpp AuxiliaryField.h FermionMatrix.cpp FermionMatrix.h HMC.cpp HMC.h CL.cpp CL.h Observables.cpp Observables.h)
add_executable(Neelix ${SOURCE_FILES} AuxiliaryField.cpp AuxiliaryField.h FermionMatrix.cpp FermionMatrix.h)
target_link_libraries(Neelix ${ARMADILLO_LIBRARIES})

set(SOURCE_FILES_UNITTST UnitTesting.cpp AuxiliaryField.cpp AuxiliaryField.h FermionMatrix.cpp FermionMatrix.h HMC.cpp HMC.h)
add_executable(unittst ${SOURCE_FILES_UNITTST})
target_link_libraries(unittst ${ARMADILLO_LIBRARIES})