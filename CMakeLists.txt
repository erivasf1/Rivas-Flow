cmake_minimum_required(VERSION 3.10)

# Project name
project(2DEulerEqSolver)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set(CMAKE_PREFIX_PATH "/apps/easybuild/software/owl-genoa/Catch2/2.13.10-GCCcore-13.3.0/lib")
# Fetch Catch2 in the system
find_package(Catch2 REQUIRED)


# Set the default build type to Debug if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# use mpicc and mpicxx to compile
SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

# Compiler flags
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-reorder")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wpedantic -Wshadow -Wconversion -fsanitize=address -lmpi")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
add_compile_definitions(OMPI_SKIP_MPICXX=1)



# Define source files
set(SOURCES
    DataManager.cpp
    EulerOperator.cpp
    ExactNozzle.cpp
    MeshAccess.hpp
    MeshGen.cpp
    Main.cpp
    Output.cpp
    TimeIntegrator.cpp
)
# Define unit testing source files
set(TESTSOURCES
    DataManager.cpp
    EulerOperator.cpp
    ExactNozzle.cpp
    MeshAccess.hpp
    MeshGen.cpp
    Output.cpp
    TimeIntegrator.cpp
    Test.cpp
)

set(OUTPUTSOURCES
    DataManager.cpp
    DiscretizationErrorOutput.cpp
    EulerOperator.cpp
    ExactNozzle.cpp
    Output.cpp
    MeshGen.cpp
    TimeIntegrator.cpp
)

# Define header files
set(HEADERS
    DataManager.h
    EulerOperator.h
    ExactNozzle.h
    MeshGen.h
    Output.h
    TimeIntegrator.h
)

set(TESTHEADERS
    DataManager.h
    EulerOperator.h
    ExactNozzle.h
    MeshGen.h
    Output.h
    TimeIntegrator.h
)

set(OUTPUTHEADERS
    DataManager.h
    EulerOperator.h
    ExactNozzle.h
    MeshGen.h
    Output.h
    TimeIntegrator.h
)

# Create the executable
add_executable(nozzle ${SOURCES} ${HEADERS})
#add_executable(testnozzle UnitTest.cpp)
add_executable(testnozzle ${TESTSOURCES} ${TESTHEADERS})
#add_executable(output ${OUTPUTSOURCES} ${OUTPUTHEADERS})

# Include directories (if needed)
#target_include_directories(NozzleSimulation PRIVATE ${CMAKE_SOURCE_DIR})

# Link libraries (if any)
#target_link_libraries(nozzle PRIVATE ${CMAKE_PREFIX_PATH})
#target_link_libraries(NozzleSimulation <LIBRARY_NAME>)
