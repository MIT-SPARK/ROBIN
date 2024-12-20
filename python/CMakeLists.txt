cmake_minimum_required(VERSION 3.18)

project(robin_python_bindings)

find_package(pybind11 CONFIG REQUIRED)

message(STATUS "Python Interpreter Version: ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}")

# fix for clang
# see: https://github.com/pybind/pybind11/issues/1818
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(robin_py PUBLIC -fsized-deallocation)
endif ()

if (DEFINED SKBUILD)
    message(STATUS "Building with Scikit-Build")
endif ()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../ ${CMAKE_CURRENT_BINARY_DIR}/robin_cpp)
else()
  message(STATUS "Performing out-of-tree build, fetching ROBIN v${CMAKE_PROJECT_VERSION} Release from Github")
  include(FetchContent)
  FetchContent_Declare(
    ext_robin PREFIX robin
    URL https://github.com/MIT-SPARK/ROBIN/archive/refs/tags/v${CMAKE_PROJECT_VERSION}.tar.gz SOURCE_SUBDIR
        cpp/patchworkpp)
  FetchContent_MakeAvailable(ext_robin)
endif()

pybind11_add_module(robin_py robin_py/robin_py.cpp)

target_link_libraries(robin_py PUBLIC robin::robin)

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  target_compile_options(robin_py PUBLIC -fsized-deallocation)
endif()

install(TARGETS robin_py DESTINATION .)
