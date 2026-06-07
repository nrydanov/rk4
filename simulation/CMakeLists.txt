cmake_minimum_required(VERSION 3.10)
project(VanDerPolSimulation)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(OpenMP REQUIRED)

set(MAIN_SIM src/demo.cpp)

set(MAIN_SYNC_ENERGY experiments/sync_energy/main.cpp)
set(MAIN_PHASE_SLOPES experiments/phase_slopes/main.cpp)
set(MAIN_MULTISTABILITY experiments/multistability/main.cpp)

add_library(vdp_core INTERFACE)

target_include_directories(vdp_core INTERFACE
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

# Sync energy experiment (collective signal power)
add_executable(vdp_sync_energy ${MAIN_SYNC_ENERGY})
target_link_libraries(vdp_sync_energy PRIVATE
    vdp_core
    yaml-cpp::yaml-cpp
    OpenMP::OpenMP_CXX
)

# Phase slopes experiment (phase difference linear regression)
add_executable(vdp_phase_slopes ${MAIN_PHASE_SLOPES})
target_link_libraries(vdp_phase_slopes PRIVATE
    vdp_core
    yaml-cpp::yaml-cpp
    OpenMP::OpenMP_CXX
)

# Multistability experiment (4 coupling types x random initial conditions)
add_executable(vdp_multistability ${MAIN_MULTISTABILITY})
target_link_libraries(vdp_multistability PRIVATE
    vdp_core
    yaml-cpp::yaml-cpp
    OpenMP::OpenMP_CXX
)

if(MSVC)
    target_compile_options(vdp_sim             PRIVATE /W4)
    target_compile_options(vdp_sync_energy     PRIVATE /W4)
    target_compile_options(vdp_phase_slopes    PRIVATE /W4)
    target_compile_options(vdp_multistability  PRIVATE /W4)
else()
    target_compile_options(vdp_sim             PRIVATE -Wall -Wextra)
    target_compile_options(vdp_sync_energy     PRIVATE -Wall -Wextra)
    target_compile_options(vdp_phase_slopes    PRIVATE -Wall -Wextra)
    target_compile_options(vdp_multistability  PRIVATE -Wall -Wextra)
    target_compile_options(vdp_sim             PRIVATE $<$<CONFIG:Release>:-O3 -march=native -DNDEBUG>)
    target_compile_options(vdp_sync_energy     PRIVATE $<$<CONFIG:Release>:-O3 -march=native -DNDEBUG>)
    target_compile_options(vdp_phase_slopes    PRIVATE $<$<CONFIG:Release>:-O3 -march=native -DNDEBUG>)
    target_compile_options(vdp_multistability  PRIVATE $<$<CONFIG:Release>:-O3 -march=native -DNDEBUG>)
endif()
