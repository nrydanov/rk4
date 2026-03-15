cmake_minimum_required(VERSION 3.10)
project(VanDerPolSimulation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(OpenMP REQUIRED)

file(GLOB_RECURSE CORE_SOURCES
    "src/solver.cpp"
    "src/forces.cpp"
)

set(MAIN_SIM src/demo.cpp)

set(MAIN_MAPS experiments/heatmap/main.cpp)

add_library(vdp_core STATIC ${CORE_SOURCES})

target_include_directories(vdp_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/external
)

include(FetchContent)

FetchContent_Declare(
  yaml-cpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
  GIT_TAG yaml-cpp-0.9.0
)
FetchContent_MakeAvailable(yaml-cpp)

# Completely YAML configured simulation
add_executable(vdp_sim ${MAIN_SIM})
target_link_libraries(vdp_sim PRIVATE
    vdp_core
    yaml-cpp::yaml-cpp
)

# Heatmap
add_executable(vdp_heatmap ${MAIN_MAPS})
target_link_libraries(vdp_heatmap PRIVATE
    vdp_core
    yaml-cpp::yaml-cpp
    OpenMP::OpenMP_CXX
)

if(MSVC)
    target_compile_options(vdp_sim  PRIVATE /W4)
    target_compile_options(vdp_heatmap PRIVATE /W4)
else()
    target_compile_options(vdp_sim  PRIVATE -Wall -Wextra)
    target_compile_options(vdp_heatmap PRIVATE -Wall -Wextra)
endif()
