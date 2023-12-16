cmake_minimum_required(VERSION 3.9)
project(RBGenerator)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/RBGeneratorCpp/src
                    ${CMAKE_CURRENT_SOURCE_DIR}/RBGeneratorCpp/ThirdParty
                    ${CMAKE_CURRENT_SOURCE_DIR}/RBGeneratorCpp/ThirdParty/pybind11
                    ${CMAKE_CURRENT_SOURCE_DIR}/RBGeneratorCpp/ThirdParty/fmt/include
                    )

set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_BINARY_DIR}/lib)

set(CMAKE_CXX_STANDARD 17)

if (UNIX)
    add_definitions("-fPIC")
endif()

find_package(OpenMP REQUIRED)

find_package(PythonInterp 3 REQUIRED)
find_package(PythonLibs 3 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

add_subdirectory(RBGeneratorCpp)
