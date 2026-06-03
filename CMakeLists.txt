cmake_minimum_required(VERSION 3.22)
set(NAME kepler-solver)
project(${NAME} LANGUAGES C)

# Build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

# Architecture string.
set(ARCH "${CMAKE_SYSTEM_PROCESSOR}" CACHE STRING "Override CPU architecture")
string(REGEX REPLACE "amd64" "x86_64" KERNEL_ARCH "${ARCH}")

# Set directory locations (allowing us to move directories easily)
set(BUILD_DIR ${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE}/${ARCH})

set(COMMON_HEADERS "./src/ks.h")
set(COMMON_SOURCES "./src/ks_bisection.c" "./src/ks_newton.c" "./src/ks_enrke.c")

# Add exe
set(EXE_NAME "ksol")
set(EXE_SOURCE "./src/main.c")
set(EXE_SOURCES ${COMMON_HEADERS} ${COMMON_SOURCES} ${EXE_SOURCE})

add_executable(${EXE_NAME} ${EXE_SOURCES})
set_target_properties(${EXE_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${BUILD_DIR}
    C_STANDARD 99
    C_STANDARD_REQUIRED ON
    C_EXTENSIONS OFF
    INTRAPROCEDURAL_OPTIMIZATION ON
)

# Target compilation options
if(MSVC)
    target_compile_options(${EXE_NAME} PRIVATE /W4 /D_CRT_SECURE_NO_WARNINGS)
    target_compile_options(${EXE_NAME} PRIVATE /fp:fast)
else()
    target_compile_options(${EXE_NAME} PRIVATE -Wall -Wextra -pedantic)
    target_compile_options(${EXE_NAME} PRIVATE -ffast-math)
endif()

# Include directories
target_include_directories(${EXE_NAME}
    PRIVATE
    ${PROJECT_SOURCE_DIR}/src
)
