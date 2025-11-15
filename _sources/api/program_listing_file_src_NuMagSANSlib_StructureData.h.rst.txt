
.. _program_listing_file_src_NuMagSANSlib_StructureData.h:

Program Listing for File NuMagSANSlib_StructureData.h
=====================================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_NuMagSANSlib_StructureData.h>` (``src/NuMagSANSlib_StructureData.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   // File         : NuMagSANSlib_MagData.h
   // Author       : Michael Philipp ADAMS, M.Sc.
   // Company      : University of Luxembourg
   // Department   : Department of Physics and Materials Sciences
   // Group        : NanoMagnetism Group
   // Group Leader : Prof. Andreas Michels
   // Version      : 23 November 2024
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
   
   struct StructureData{
   
      float *RotMat;   // rotation matrix R = [R0, R3, R6; R1, R4, R7; R2, R5, R8]
      float*x;     // position data x
      float*y;     // position data y
      float*z;     // position data z
      unsigned long int* K;    // number of objects
   
   };
   
   void allocate_StructureDataRAM(StructureData* StructData, \
                                  StructDataProperties* StructDataProp, \
                                  InputFileData* InputData){
   
   
        unsigned long int K = StructDataProp->Number_Of_Elements;
   
        StructData->x = (float*) malloc(K*sizeof(float));
        StructData->y = (float*) malloc(K*sizeof(float));
        StructData->z = (float*) malloc(K*sizeof(float));
        StructData->K = (unsigned long int*) malloc(sizeof(unsigned long int));
        StructData->RotMat = (float*) malloc(9*sizeof(float));
   
        if (!StructData->x || !StructData->y || !StructData->z ||
            !StructData->K || !StructData->RotMat) {
                perror("Memory allocation failed");
                exit(EXIT_FAILURE);
        }
   
        *StructData->K = K;
   
        for(unsigned long int i = 0; i < K; i++){
           StructData->x[i] = 0.0;
           StructData->y[i] = 0.0;
           StructData->z[i] = 0.0;
        }  
   }
   
   
   void allocate_StructureDataGPU(StructureData* StructData, \
                               StructureData* StructData_gpu){
   
        unsigned long int K = *StructData->K;
       
        cudaMalloc(&StructData_gpu->x, K*sizeof(float));
        cudaMalloc(&StructData_gpu->y, K*sizeof(float));
        cudaMalloc(&StructData_gpu->z, K*sizeof(float));
        cudaMalloc(&StructData_gpu->K, sizeof(unsigned long int));
        cudaMalloc(&StructData_gpu->RotMat, 9*sizeof(float));
   
        cudaMemcpy(StructData_gpu->K, StructData->K, sizeof(unsigned long int), cudaMemcpyHostToDevice);
        cudaMemcpy(StructData_gpu->RotMat, StructData->RotMat, 9*sizeof(float), cudaMemcpyHostToDevice);
           
       LogSystem::write("");
            // copy data from Host to Device
       LogSystem::write("copy data from RAM to GPU...");
        cudaMemcpy(StructData_gpu->x, StructData->x, K*sizeof(float), cudaMemcpyHostToDevice);
       LogSystem::write("   x done...");
        cudaMemcpy(StructData_gpu->y, StructData->y, K*sizeof(float), cudaMemcpyHostToDevice);
       LogSystem::write("   y done...");
        cudaMemcpy(StructData_gpu->z, StructData->z, K*sizeof(float), cudaMemcpyHostToDevice);
       LogSystem::write("   z done...");
       LogSystem::write("");
       LogSystem::write("data transfer finished...");
       LogSystem::write("");
   
       cudaMalloc(&StructData_gpu, sizeof(StructureData));
       cudaMemcpy(StructData_gpu, StructData, sizeof(StructureData), cudaMemcpyHostToDevice);
        
   }
   
   
   void read_StructureData(StructureData* StructData, \
                           StructDataProperties* StructDataProp, \
                           InputFileData* InputData){
   
       LogSystem::write("");
       LogSystem::write("read StructData...");
   
        //unsigned long int K = *StructData->K;
   
        for(int i = 0; i < 9; i++){
           StructData->RotMat[i] = InputData->RotMat[i];
        }
   
        string filename;
        unsigned long int n = 0;
        float x_buf, y_buf, z_buf;
        ifstream fin;
   
   
        filename = StructDataProp->GlobalFilePath;
       LogSystem::write(filename);
   
   
        fin.open(filename);
        n = 0;
        // read in the data
        while(fin >> x_buf >> y_buf >> z_buf){
           StructData->x[n] = x_buf * InputData->XYZ_Unit_Factor;
           StructData->y[n] = y_buf * InputData->XYZ_Unit_Factor;
           StructData->z[n] = z_buf * InputData->XYZ_Unit_Factor;
           n += 1;
       }
   
       fin.close();
       LogSystem::write("read (x,y,z) StructData finished...");
       
   }
   
   void init_StructureData(StructureData* StructData, \
                           StructureData* StructData_gpu, \
                           StructDataProperties* StructDataProp, \
                           InputFileData* InputData){
   
       allocate_StructureDataRAM(StructData, StructDataProp, InputData);
       read_StructureData(StructData, StructDataProp, InputData);
       allocate_StructureDataGPU(StructData, StructData_gpu);
      
   }
   
   void free_StructureData(StructureData *StructData, \
                           StructureData *StructData_gpu){
   
       LogSystem::write("free StructData...");
        free(StructData->x);
        free(StructData->y);
        free(StructData->z);
        free(StructData->K);
        free(StructData->RotMat);
   
        cudaDeviceSynchronize();
        cudaFree(StructData_gpu->x);
        cudaDeviceSynchronize();
        cudaFree(StructData_gpu->y);
        cudaDeviceSynchronize();
        cudaFree(StructData_gpu->z);
        cudaDeviceSynchronize();
        cudaFree(StructData_gpu->K);
        cudaDeviceSynchronize();
        cudaFree(StructData_gpu->RotMat);
        cudaDeviceSynchronize();
   
       LogSystem::write("free StructData finished.");
   
   }
   
   void disp_StructureData(StructureData *StructData){
       for(int k=0; k < *StructData->K; k++){
           cout << StructData->x[k] << " "\
                << StructData->y[k] << " "\
                << StructData->z[k] << "\n";
       }
   }
