cmake_minimum_required(VERSION 3.10)
project(NuMagSANS)

enable_language(CUDA)
set (CUDA_STANDARD 11)

add_executable(NuMagSANS src/NuMagSANS.cu)
