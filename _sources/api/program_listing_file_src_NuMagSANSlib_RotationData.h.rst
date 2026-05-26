
.. _program_listing_file_src_NuMagSANSlib_RotationData.h:

Program Listing for File NuMagSANSlib_RotationData.h
====================================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_NuMagSANSlib_RotationData.h>` (``src/NuMagSANSlib_RotationData.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

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
   
   
   // The RotationData structure is aimed to allow rotations of individual objects 
   // in a ensemble of objects. This means the data handled through this structure
   // are diffent from the rotation angles in the InputDataFile, which are aimed for 
   // a global rotation of the complete system, where the angles here allow object-wise rotation.
   
   struct RotationData{
   
       // The rotation angles correspond to a zyz-rotation matrix combination
       // The convention of the three rotation angles is defined as follows:
       // R(alpha, beta, gamma) = R_z(alpha) * Ry(beta) * R_z(gamma)
       // Rotations are assumed to be performed as M = R * m, X = R * x
       // the angles are assumed in radian format
       float* alpha;  
       float* beta;
       float* gamma;
   
       // The RotMat variable compacts the three rotation angles to a direct
       // rotation matrix variable with the indice convention:
       //R = [R0, R3, R6; 
       //     R1, R4, R7; 
       //     R2, R5, R8]
       float* RotMat; 
   
       unsigned long int* K;   // number of objects
   
   };
   
   
   void allocate_RotationDataRAM(RotationData* RotData, \
                                  RotDataProperties* RotDataProp, \
                                  InputFileData* InputData){
   
   
        unsigned long int K = RotDataProp->Number_Of_Elements;
   
        RotData->alpha = (float*) malloc(K*sizeof(float));
        RotData->beta = (float*) malloc(K*sizeof(float));
        RotData->gamma = (float*) malloc(K*sizeof(float));
        RotData->K = (unsigned long int*) malloc(sizeof(unsigned long int));
        RotData->RotMat = (float*) malloc(9*K*sizeof(float));
   
        if (!RotData->alpha || !RotData->beta || !RotData->gamma ||
            !RotData->K || !RotData->RotMat) {
                perror("Memory allocation failed");
                exit(EXIT_FAILURE);
        }
   
        *RotData->K = K;
   
        for(unsigned long int i = 0; i < K; i++){
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
           RotData->RotMat[9*i+0] = 1; 
           RotData->RotMat[9*i+1] = 0;
           RotData->RotMat[9*i+2] = 0;
   
           RotData->RotMat[9*i+3] = 0;
           RotData->RotMat[9*i+4] = 1;
           RotData->RotMat[9*i+5] = 0;
   
           RotData->RotMat[9*i+6] = 0;
           RotData->RotMat[9*i+7] = 0;
           RotData->RotMat[9*i+8] = 1;
        }  
   }
   
   void allocate_RotationDataGPU(RotationData* RotData, \
                                 RotationData* RotData_gpu){
   
        unsigned long int K = *RotData->K;
       
       cudaMalloc(&RotData_gpu->alpha, K*sizeof(float));
       cudaMalloc(&RotData_gpu->beta, K*sizeof(float));
       cudaMalloc(&RotData_gpu->gamma, K*sizeof(float));
       cudaMalloc(&RotData_gpu->K, sizeof(unsigned long int));
       cudaMalloc(&RotData_gpu->RotMat, 9*K*sizeof(float));
   
       cudaMemcpy(RotData_gpu->K, RotData->K, sizeof(unsigned long int), cudaMemcpyHostToDevice);
       cudaMemcpy(RotData_gpu->RotMat, RotData->RotMat, 9*K*sizeof(float), cudaMemcpyHostToDevice);
           
       LogSystem::write("");
            // copy data from Host to Device
       LogSystem::write("copy data from RAM to GPU...");
        cudaMemcpy(RotData_gpu->alpha, RotData->alpha, K*sizeof(float), cudaMemcpyHostToDevice);
       LogSystem::write("   alpha done...");
        cudaMemcpy(RotData_gpu->beta, RotData->beta, K*sizeof(float), cudaMemcpyHostToDevice);
       LogSystem::write("   beta done...");
        cudaMemcpy(RotData_gpu->gamma, RotData->gamma, K*sizeof(float), cudaMemcpyHostToDevice);
       LogSystem::write("   gamma done...");
       LogSystem::write("");
       LogSystem::write("data transfer finished...");
       LogSystem::write("");
   
       cudaMalloc(&RotData_gpu, sizeof(RotationData));
       cudaMemcpy(RotData_gpu, RotData, sizeof(RotationData), cudaMemcpyHostToDevice);
        
   }
   
   
   // the followin section defins rotation matrix operations 
   // later these operations can be separated to a new header file
   // to make handle the code more separated
   void RotationMatrix_z(float alpha, float* Rz){
   
       // This function generates a z-Rotation matrix based on the angle
       // alpha as a single column array. The ordering of the single column
       // array is understood as follows: 
       // Rz = [R0, R3, R6
       //       R1, R4, R7, 
       //       R2, R5, R8] 
   
       float c = cos(alpha);
       float s = sin(alpha);
   
       Rz[0] = c;
       Rz[1] = s;
       Rz[2] = 0.0;
   
       Rz[3] = -s;
       Rz[4] = c;
       Rz[5] = 0.0;
   
       Rz[6] = 0.0;
       Rz[7] = 0.0;
       Rz[8] = 1.0;
    
   }
   
   void RotationMatrix_y(float alpha, float* Ry){
   
       // This function generates a y-Rotation matrix based on the angle
       // alpha as a single column array. The ordering of the single column
       // array is understood as follows: 
       // Ry = [R0, R3, R6
       //       R1, R4, R7, 
       //       R2, R5, R8] 
   
       float c = cos(alpha);
       float s = sin(alpha);
   
       Ry[0] = c;
       Ry[1] = 0.0;
       Ry[2] = -s;
   
       Ry[3] = 0.0;
       Ry[4] = 1.0;
       Ry[5] = 0.0;
   
       Ry[6] = s;
       Ry[7] = 0.0;
       Ry[8] = c;
    
   }
   
   void LeftMultiply3x3(float* R1, const float* R2) {
       // Updates R1 in place:
       //
       //     R1 <- R2 * R1
       //
       // Matrix storage is column-major:
       //
       //     R(row, col) = R[row + 3 * col]
       //
       // Therefore:
       //
       //     R1_new(row, col) = sum_k R2(row, k) * R1_old(k, col)
   
       float T[9];
   
       for (int col = 0; col < 3; ++col) {
           for (int row = 0; row < 3; ++row) {
               T[row + 3 * col] =
                   R2[row + 3 * 0] * R1[0 + 3 * col] +
                   R2[row + 3 * 1] * R1[1 + 3 * col] +
                   R2[row + 3 * 2] * R1[2 + 3 * col];
           }
       }
   
       for (int i = 0; i < 9; ++i) {
           R1[i] = T[i];
       }
   }
   
   void Multiply_RotmatZYZ_3x3(float alpha, float beta, float gamma, float* RotMat){
   
       // Assumption:
       // RotMat is initialized as identity matrix before entering this function.
       //
       // Iterative update:
       // RotMat <- Rz(gamma) * RotMat
       // RotMat <- Ry(beta)  * RotMat
       // RotMat <- Rz(alpha) * RotMat
       //
       // Final result:
       // RotMat = Rz(alpha) * Ry(beta) * Rz(gamma)
   
       float Rz1[9];
       float Ry1[9];
       float Rz2[9];
   
       RotationMatrix_z(alpha, Rz1);
       RotationMatrix_y(beta, Ry1);
       RotationMatrix_z(gamma, Rz2);
   
       LeftMultiply3x3(RotMat, Rz2);
       LeftMultiply3x3(RotMat, Ry1);
       LeftMultiply3x3(RotMat, Rz1);
   }
   
   void RotMat_select(unsigned long int idx, float* RotMat, float*RotMat_idx){
   
       // this function extracts a single rotation matrix from the large
       // rotation matrix vector at given index idx
   
       RotMat_idx[0] = RotMat[9*idx+0];
       RotMat_idx[1] = RotMat[9*idx+1];
       RotMat_idx[2] = RotMat[9*idx+2];
       RotMat_idx[3] = RotMat[9*idx+3];
       RotMat_idx[4] = RotMat[9*idx+4];
       RotMat_idx[5] = RotMat[9*idx+5];
       RotMat_idx[6] = RotMat[9*idx+6];
       RotMat_idx[7] = RotMat[9*idx+7];
       RotMat_idx[8] = RotMat[9*idx+8];
   
   }
   
   void RotMat_store(unsigned long int idx, float* RotMat, float*RotMat_idx){
   
       // this function stores a single rotation matrix to the large
       // rotation matrix vector at given index idx
   
       RotMat[9*idx+0] = RotMat_idx[0];
       RotMat[9*idx+1] = RotMat_idx[1];
       RotMat[9*idx+2] = RotMat_idx[2];
       RotMat[9*idx+3] = RotMat_idx[3];
       RotMat[9*idx+4] = RotMat_idx[4];
       RotMat[9*idx+5] = RotMat_idx[5];
       RotMat[9*idx+6] = RotMat_idx[6];
       RotMat[9*idx+7] = RotMat_idx[7];
       RotMat[9*idx+8] = RotMat_idx[8];
       
   }
   
   
   // here the rotation matrix operations section ends
   
   
   void read_RotationData(RotationData* RotData, \
                          RotDataProperties* RotDataProp, \
                          InputFileData* InputData){
   
       LogSystem::write("");
       LogSystem::write("read RotationData...");
   
       unsigned long int K = *RotData->K;
   
       string filename;
       unsigned long int n = 0;
       float alpha_buf, beta_buf, gamma_buf;
       ifstream fin;
   
       filename = RotDataProp->GlobalFilePath;
       LogSystem::write(filename);
   
   
        fin.open(filename);
        n = 0;
        // read in the data
        while(fin >> alpha_buf >> beta_buf >> gamma_buf){
           RotData->alpha[n] = alpha_buf;
           RotData->beta[n] = beta_buf;
           RotData->gamma[n] = gamma_buf;
           n += 1;
       }
   
       float RotMat_buf[9];
       for(unsigned long int k=0; k<K; k++){
   
           RotMat_select(k, RotData->RotMat, RotMat_buf);
   
           Multiply_RotmatZYZ_3x3(RotData->alpha[k], 
                                  RotData->beta[k], 
                                  RotData->gamma[k], 
                                  RotMat_buf);
   
           RotMat_store(k, RotData->RotMat, RotMat_buf);
   
       }
   
       fin.close();
       LogSystem::write("read (alpha, beta, gamma) RotationData finished...");
       
   }
   
   void init_RotationData(RotationData* RotData, \
                           RotationData* RotData_gpu, \
                           RotDataProperties* RotDataProp, \
                           InputFileData* InputData){
   
       allocate_RotationDataRAM(RotData, RotDataProp, InputData);
       read_RotationData(RotData, RotDataProp, InputData);
       allocate_RotationDataGPU(RotData, RotData_gpu);
      
   }
   
   void free_RotationData(RotationData *RotData,
                          RotationData *RotData_gpu){
   
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
   
