cmake_minimum_required(VERSION 3.20)
project(lab1 C)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include_directories(include)

file(GLOB SOURCES "src/*.c")

add_executable(main ${SOURCES})

target_link_libraries(main m)
