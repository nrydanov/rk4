cmake_minimum_required(VERSION 3.10)
project(VanDerPolSimulation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(yaml-cpp REQUIRED)

file(GLOB_RECURSE SOURCES "src/*.cpp")

add_executable(vdp_sim ${SOURCES})

target_include_directories(vdp_sim PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/external
)
target_link_libraries(vdp_sim PRIVATE yaml-cpp::yaml-cpp)

if(MSVC)
    target_compile_options(vdp_sim PRIVATE /W4)
else()
    target_compile_options(vdp_sim PRIVATE -Wall -Wextra)
endif()
