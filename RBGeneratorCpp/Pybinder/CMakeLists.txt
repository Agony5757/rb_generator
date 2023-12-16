cmake_minimum_required(VERSION 3.1)
project(RBGeneratorCpp)

# set(LIB_RBGenerator ${PROJECT_NAME})
file(GLOB_RECURSE RBGeneratorWrapper_SRC_CPP *.cpp)

pybind11_add_module(${PROJECT_NAME} ${RBGeneratorWrapper_SRC_CPP})

target_link_libraries(${PROJECT_NAME} PRIVATE fmt)

target_compile_definitions(${PROJECT_NAME} PRIVATE VERSION_INFO=${RBGenerator_VERSION_INFO})