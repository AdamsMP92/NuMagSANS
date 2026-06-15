cmake_minimum_required(VERSION 3.10)
project(NuMagSANS LANGUAGES CXX)

option(NUMAGSANS_ENABLE_HDF5 "Enable optional HDF5 output backend" OFF)

# -----------------------------------------------------------
# Optional CUDA support (detect automatically)
# -----------------------------------------------------------
include(CheckLanguage)
check_language(CUDA)

if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    message(STATUS "✅ CUDA compiler found: ${CMAKE_CUDA_COMPILER}")

    # Set CUDA standard
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_EXTENSIONS OFF)

    # Make sure host compiler also uses C++17
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)

    # Add executable
    add_executable(NuMagSANS src/NuMagSANS.cu)

    if(NUMAGSANS_ENABLE_HDF5)
        find_package(HDF5 REQUIRED COMPONENTS C)
        target_include_directories(NuMagSANS PRIVATE ${HDF5_INCLUDE_DIRS})
        target_link_libraries(NuMagSANS PRIVATE ${HDF5_LIBRARIES})
        target_compile_definitions(NuMagSANS PRIVATE NUMAGSANS_ENABLE_HDF5)
        message(STATUS "✅ HDF5 output backend enabled")
    else()
        message(STATUS "HDF5 output backend disabled; CSV output remains the default")
    endif()

    # (optional) Additional flags or architecture targets
    # set_target_properties(NuMagSANS PROPERTIES CUDA_ARCHITECTURES "60;70;75;80")

else()
    message(WARNING "⚠️ CUDA compiler not found — skipping CUDA build.")
endif()

# -----------------------------------------------------------
# General build settings
# -----------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
