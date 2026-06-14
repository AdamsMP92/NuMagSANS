#pragma once

inline bool CheckCUDA(cudaError_t err, const std::string& context) {
    if (err != cudaSuccess) {
        LogSystem::write("CUDA error in " + context + ": " + cudaGetErrorString(err));
        return false;
    }

    return true;
}

inline bool CheckCUDALastError(const std::string& context) {
    bool success = true;

    success = CheckCUDA(cudaGetLastError(), context) && success;

    success = CheckCUDA(cudaDeviceSynchronize(), context + " synchronization") && success;

    return success;
}

inline bool CheckCUDAKernelRun(const std::string& kernel_name) {
    bool success = true;

    success = CheckCUDA(cudaGetLastError(), "kernel launch " + kernel_name) && success;

    success = CheckCUDA(cudaDeviceSynchronize(), "kernel execution " + kernel_name) && success;

    return success;
}
