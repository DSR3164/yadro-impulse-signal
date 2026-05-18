cmake_minimum_required(VERSION 3.24)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    if(NOT MSVC)
        message(STATUS "Set debug mode")
        set(SANITY_FLAGS -fsanitize=address,undefined -fno-omit-frame-pointer -g)
        add_compile_options(${SANITY_FLAGS})
        add_link_options(${SANITY_FLAGS})
    endif()
endif()

project(YadroTest LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_library(core STATIC
    ${PROJECT_SOURCE_DIR}/src/dsp.cpp
)

target_include_directories(core
    PUBLIC ${PROJECT_SOURCE_DIR}/include
)

target_compile_options(core PRIVATE -mavx2 -ffast-math -Wall)

add_executable(main main.cpp)

target_link_libraries(main PRIVATE core)
