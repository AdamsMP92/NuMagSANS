cmake_minimum_required(VERSION 3.10)
project(NuMagSANS LANGUAGES CXX)

# -----------------------------------------------------------
# Optional CUDA support (detect automatically)
# -----------------------------------------------------------
include(CheckLanguage)
check_language(CUDA)

if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    message(STATUS "✅ CUDA compiler found: ${CMAKE_CUDA_COMPILER}")

    # Set CUDA standard
    set(CMAKE_CUDA_STANDARD 11)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_EXTENSIONS OFF)

    # Add executable
    add_executable(NuMagSANS src/NuMagSANS.cu)

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
