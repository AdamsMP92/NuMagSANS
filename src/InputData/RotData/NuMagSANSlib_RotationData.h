// File         : NuMagSANSlib_RotationData.h
// Author       : Michael Philipp ADAMS, M.Sc.
// Company      : University of Luxembourg
// Department   : Department of Physics and Materials Sciences
// Group        : NanoMagnetism Group
// Group Leader : Prof. Andreas Michels
// Version      : 25 May 2026
// OS           : Linux Ubuntu
// Language     : CUDA C++

#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <string>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <stdexcept>
#include <math.h>
#include <chrono>
#include <dirent.h>
#include <unistd.h>

using namespace std;

// The RotationData structure allows rotations of individual objects in an
// ensemble. These rotations are distinct from the global two-angle rotation in
// the input configuration.

struct RotationData {

    // The three angles follow InputData::RotDataConvention. The default is the
    // historical Z-Y-Z convention
    // R(alpha, beta, gamma) = R_z(alpha) * R_y(beta) * R_z(gamma).
    // Rotations are active and act on column vectors as M = R * m and X = R * x.
    // Angles are expressed in radians.
    float* alpha;
    float* beta;
    float* gamma;

    // The RotMat variable compacts the three rotation angles to a direct
    // rotation matrix variable with the indice convention:
    // R = [R0, R3, R6;
    //     R1, R4, R7;
    //     R2, R5, R8]
    float* RotMat;

    unsigned long int* K; // number of objects
};

void allocate_RotationDataRAM(RotationData* RotData, RotDataProperties* RotDataProp, InputFileData* InputData) {

    unsigned long int K = RotDataProp->Number_Of_Elements;

    RotData->alpha = (float*)malloc(K * sizeof(float));
    RotData->beta = (float*)malloc(K * sizeof(float));
    RotData->gamma = (float*)malloc(K * sizeof(float));
    RotData->K = (unsigned long int*)malloc(sizeof(unsigned long int));
    RotData->RotMat = (float*)malloc(9 * K * sizeof(float));

    if (!RotData->alpha || !RotData->beta || !RotData->gamma || !RotData->K || !RotData->RotMat) {
        perror("Memory allocation failed");
        exit(EXIT_FAILURE);
    }
}

void DefaultSet_RotationDataRAM(RotationData* RotData, RotDataProperties* RotDataProp, InputFileData* InputData) {

    unsigned long int K = RotDataProp->Number_Of_Elements;

    *RotData->K = K;

    for (unsigned long int i = 0; i < K; i++) {
        // initialize angles with zeros
        RotData->alpha[i] = 0.0;
        RotData->beta[i] = 0.0;
        RotData->gamma[i] = 0.0;

        // initialize rotation matrices as identity matrices
        // for one set of angles alpha, beta, gamma at index i
        // the rotation matrix is stored as a single column array
        // with 9 entries following the index logic
        // R = [R_{0,i}, R_{3,i}, R_{6,i},
        //      R_{1,i}, R_{4,i}, R_{7,i},
        //      R_{2,i}, R_{5,i}, R_{8,i}]
        //   = [1, 0, 0,
        //      0, 1, 0,
        //      0, 0, 1]
        RotData->RotMat[9 * i + 0] = 1;
        RotData->RotMat[9 * i + 1] = 0;
        RotData->RotMat[9 * i + 2] = 0;

        RotData->RotMat[9 * i + 3] = 0;
        RotData->RotMat[9 * i + 4] = 1;
        RotData->RotMat[9 * i + 5] = 0;

        RotData->RotMat[9 * i + 6] = 0;
        RotData->RotMat[9 * i + 7] = 0;
        RotData->RotMat[9 * i + 8] = 1;
    }
}

void allocate_RotationDataGPU(RotationData* RotData, RotationData* RotData_gpu) {

    unsigned long int K = *RotData->K;

    cudaMalloc(&RotData_gpu->alpha, K * sizeof(float));
    cudaMalloc(&RotData_gpu->beta, K * sizeof(float));
    cudaMalloc(&RotData_gpu->gamma, K * sizeof(float));
    cudaMalloc(&RotData_gpu->K, sizeof(unsigned long int));
    cudaMalloc(&RotData_gpu->RotMat, 9 * K * sizeof(float));

    cudaMemcpy(RotData_gpu->K, RotData->K, sizeof(unsigned long int), cudaMemcpyHostToDevice);
    cudaMemcpy(RotData_gpu->RotMat, RotData->RotMat, 9 * K * sizeof(float), cudaMemcpyHostToDevice);

    LogSystem::write("");
    LogSystem::write("copy data from RAM to GPU...");
    cudaMemcpy(RotData_gpu->alpha, RotData->alpha, K * sizeof(float), cudaMemcpyHostToDevice);
    LogSystem::write("   angle_1 done...");
    cudaMemcpy(RotData_gpu->beta, RotData->beta, K * sizeof(float), cudaMemcpyHostToDevice);
    LogSystem::write("   angle_2 done...");
    cudaMemcpy(RotData_gpu->gamma, RotData->gamma, K * sizeof(float), cudaMemcpyHostToDevice);
    LogSystem::write("   angle_3 done...");
    LogSystem::write("");
    LogSystem::write("data transfer finished...");
    LogSystem::write("");

    cudaMalloc(&RotData_gpu, sizeof(RotationData));
    cudaMemcpy(RotData_gpu, RotData, sizeof(RotationData), cudaMemcpyHostToDevice);
}

void copy_RotationDataRAM2GPU(RotationData* RotData, RotationData* RotData_gpu) {

    unsigned long int K = *RotData->K;

    cudaMemcpy(RotData_gpu->K, RotData->K, sizeof(unsigned long int), cudaMemcpyHostToDevice);
    cudaMemcpy(RotData_gpu->RotMat, RotData->RotMat, 9 * K * sizeof(float), cudaMemcpyHostToDevice);

    LogSystem::write("");
    LogSystem::write("copy data from RAM to GPU...");
    cudaMemcpy(RotData_gpu->alpha, RotData->alpha, K * sizeof(float), cudaMemcpyHostToDevice);
    LogSystem::write("   angle_1 done...");
    cudaMemcpy(RotData_gpu->beta, RotData->beta, K * sizeof(float), cudaMemcpyHostToDevice);
    LogSystem::write("   angle_2 done...");
    cudaMemcpy(RotData_gpu->gamma, RotData->gamma, K * sizeof(float), cudaMemcpyHostToDevice);
    LogSystem::write("   angle_3 done...");
    LogSystem::write("");
    LogSystem::write("data transfer finished...");
    LogSystem::write("");
}

void RotMat_select(unsigned long int idx, float* RotMat, float* RotMat_idx) {

    // this function extracts a single rotation matrix from the large
    // rotation matrix vector at given index idx

    RotMat_idx[0] = RotMat[9 * idx + 0];
    RotMat_idx[1] = RotMat[9 * idx + 1];
    RotMat_idx[2] = RotMat[9 * idx + 2];
    RotMat_idx[3] = RotMat[9 * idx + 3];
    RotMat_idx[4] = RotMat[9 * idx + 4];
    RotMat_idx[5] = RotMat[9 * idx + 5];
    RotMat_idx[6] = RotMat[9 * idx + 6];
    RotMat_idx[7] = RotMat[9 * idx + 7];
    RotMat_idx[8] = RotMat[9 * idx + 8];
}

void RotMat_store(unsigned long int idx, float* RotMat, float* RotMat_idx) {

    // this function stores a single rotation matrix to the large
    // rotation matrix vector at given index idx

    RotMat[9 * idx + 0] = RotMat_idx[0];
    RotMat[9 * idx + 1] = RotMat_idx[1];
    RotMat[9 * idx + 2] = RotMat_idx[2];
    RotMat[9 * idx + 3] = RotMat_idx[3];
    RotMat[9 * idx + 4] = RotMat_idx[4];
    RotMat[9 * idx + 5] = RotMat_idx[5];
    RotMat[9 * idx + 6] = RotMat_idx[6];
    RotMat[9 * idx + 7] = RotMat_idx[7];
    RotMat[9 * idx + 8] = RotMat_idx[8];
}

void read_RotationData(RotationData* RotData, RotDataProperties* RotDataProp, InputFileData* InputData) {

    LogSystem::write("");
    LogSystem::write("read RotationData...");

    unsigned long int K = *RotData->K;

    string filename;
    unsigned long int n = 0;
    float angle_1_buf, angle_2_buf, angle_3_buf;
    ifstream fin;

    filename = RotDataProp->GlobalFilePath;
    LogSystem::write(filename);

    fin.open(filename);
    n = 0;
    // read in the data
    while (fin >> angle_1_buf >> angle_2_buf >> angle_3_buf) {
        RotData->alpha[n] = angle_1_buf;
        RotData->beta[n] = angle_2_buf;
        RotData->gamma[n] = angle_3_buf;
        n += 1;
    }

    float RotMat_buf[9];
    for (unsigned long int k = 0; k < K; k++) {

        RotMat_select(k, RotData->RotMat, RotMat_buf);

        Multiply_RotmatEuler_3x3(RotData->alpha[k], RotData->beta[k], RotData->gamma[k], InputData->RotDataConvention,
                                 RotMat_buf);

        RotMat_store(k, RotData->RotMat, RotMat_buf);
    }

    fin.close();
    LogSystem::write("read (angle_1, angle_2, angle_3) RotationData finished...");
}

void init_RotationData(RotationData* RotData, RotationData* RotData_gpu, RotDataProperties* RotDataProp,
                       InputFileData* InputData) {

    allocate_RotationDataRAM(RotData, RotDataProp, InputData);
    DefaultSet_RotationDataRAM(RotData, RotDataProp, InputData);
    read_RotationData(RotData, RotDataProp, InputData);
    allocate_RotationDataGPU(RotData, RotData_gpu);
}

void init_RotationDataMemory(RotationData* RotData, RotationData* RotData_gpu, RotDataProperties* RotDataProp,
                             InputFileData* InputData) {

    allocate_RotationDataRAM(RotData, RotDataProp, InputData);
    DefaultSet_RotationDataRAM(RotData, RotDataProp, InputData);
    allocate_RotationDataGPU(RotData, RotData_gpu);
}

void new_read_RotationData(RotationData* RotData, RotationData* RotData_gpu, RotDataProperties* RotDataProp,
                           InputFileData* InputData) {

    DefaultSet_RotationDataRAM(RotData, RotDataProp, InputData);
    read_RotationData(RotData, RotDataProp, InputData);
    copy_RotationDataRAM2GPU(RotData, RotData_gpu);
}

void free_RotationData(RotationData* RotData, RotationData* RotData_gpu) {

    LogSystem::write("free RotationData...");

    free(RotData->alpha);
    free(RotData->beta);
    free(RotData->gamma);
    free(RotData->K);
    free(RotData->RotMat);

    cudaDeviceSynchronize();
    cudaFree(RotData_gpu->alpha);

    cudaDeviceSynchronize();
    cudaFree(RotData_gpu->beta);

    cudaDeviceSynchronize();
    cudaFree(RotData_gpu->gamma);

    cudaDeviceSynchronize();
    cudaFree(RotData_gpu->K);

    cudaDeviceSynchronize();
    cudaFree(RotData_gpu->RotMat);

    cudaDeviceSynchronize();

    LogSystem::write("free RotationData finished.");
}
