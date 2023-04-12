cmake_minimum_required(VERSION 3.10)

project(robin_python_bindings)

pybind11_add_module(robin_py robin_py/robin_py.cpp)

message(STATUS "Python Interpreter Version: ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}")

target_link_libraries(robin_py PUBLIC robin)

# fix for clang
# see: https://github.com/pybind/pybind11/issues/1818
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(robin_py PUBLIC -fsized-deallocation)
endif ()

# make sure to output the build file to robin_py folder
SET_TARGET_PROPERTIES(robin_py
        PROPERTIES
        PREFIX ""
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/robin_py"
        )

# copy package __init__.py file
configure_file(robin_py/__init__.py
        ${CMAKE_CURRENT_BINARY_DIR}/robin_py/__init__.py
        )

# copy setup.py file
configure_file(setup.py.in
        ${CMAKE_CURRENT_BINARY_DIR}/setup.py
        )

file(COPY .
        DESTINATION .
        FILES_MATCHING
        PATTERN *.py)