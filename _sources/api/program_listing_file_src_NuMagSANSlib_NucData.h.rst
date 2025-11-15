
.. _program_listing_file_src_NuMagSANSlib_NucData.h:

Program Listing for File NuMagSANSlib_NucData.h
===============================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_NuMagSANSlib_NucData.h>` (``src/NuMagSANSlib_NucData.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   // File         : NuMagSANSlib_NucData.h
   // Author       : Michael Philipp ADAMS, M.Sc.
   // Company      : University of Luxembourg
   // Department   : Department of Physics and Materials Sciences
   // Group        : NanoMagnetism Group
   // Group Leader : Prof. Andreas Michels
   // Version      : 26 November 2024
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
   
   struct NuclearData{
   
       float *RotMat;   // rotation matrix R = [R0, R3, R6; R1, R4, R7; R2, R5, R8]
       float*x;     // position data x
       float*y;     // position data y
       float*z;     // position data z
       float*Nuc;    // nuclear density data
   
       unsigned long int* K;    // number of orientations
       unsigned long int* TotalAtomNumber;  // total atom number
   
       unsigned int* NumberOfElements;
       unsigned long int* N_cum;   // number of atoms in k-th object
       unsigned long int* N_act;   // cummulated number of atoms from object to object
       unsigned long int* N_avg;   // number of nuclear points
   
   };
   
   void allocate_NuclearDataRAM(NuclearData* NucData, \
                                NucDataProperties* NucDataProp, \
                                int NucData_File_Index){
       LogSystem::write("data file index: " + std::to_string(NucData_File_Index));
       LogSystem::write("total atom number: " + std::to_string(NucDataProp->TotalAtomNumber[0]));
   
   
   
            unsigned long int K = NucDataProp->Number_Of_SubFolders;
            unsigned long int TotalAtomNumber = NucDataProp->TotalAtomNumber[NucData_File_Index-1];
      
            NucData->x = (float*) malloc(TotalAtomNumber*sizeof(float));
            NucData->y = (float*) malloc(TotalAtomNumber*sizeof(float));
            NucData->z = (float*) malloc(TotalAtomNumber*sizeof(float));
            NucData->Nuc = (float*) malloc(TotalAtomNumber*sizeof(float));
            NucData->K = (unsigned long int*) malloc(sizeof(unsigned long int));
            NucData->TotalAtomNumber = (unsigned long int*) malloc(sizeof(unsigned long int));
            NucData->RotMat = (float*) malloc(9*sizeof(float));
            NucData->NumberOfElements = (unsigned int*) malloc(K*sizeof(int));
            NucData->N_cum = (unsigned long int*) malloc(K*sizeof(unsigned long int));
            NucData->N_act = (unsigned long int*) malloc(K*sizeof(unsigned long int));
            NucData->N_avg = (unsigned long int*) malloc(sizeof(unsigned long int));
   
           if(!NucData->x || !NucData->y || !NucData->z || !NucData->Nuc ||
              !NucData->N_avg || !NucData->K || !NucData->TotalAtomNumber || !NucData->RotMat){
                   perror("Memory alloation failed");
                   exit(EXIT_FAILURE);
           }
      
            *NucData->N_avg = 0;
            *NucData->K = K;
            *NucData->TotalAtomNumber = TotalAtomNumber;
   
            for(unsigned long int i = 0; i < TotalAtomNumber; i++){
               NucData->x[i] = 0.0;
               NucData->y[i] = 0.0;
               NucData->z[i] = 0.0;
               NucData->Nuc[i] = 0.0;
            }
   
            for(unsigned int i = 0; i < K; i++){
                NucData->NumberOfElements[i] = NucDataProp->NumberOfElements[i][NucData_File_Index-1];
                NucData->N_cum[i] = 0;
                NucData->N_act[i] = 0;
           }
   
   }
   
   
    void allocate_NuclearDataGPU(NuclearData* NucData, \
                                 NuclearData* NucData_gpu,\
                                 int NucData_File_Index){
     
           unsigned long int TotalAtomNumber = *NucData->TotalAtomNumber;
     
           cudaMalloc(&NucData_gpu->x, TotalAtomNumber*sizeof(float));
           cudaMalloc(&NucData_gpu->y, TotalAtomNumber*sizeof(float));
           cudaMalloc(&NucData_gpu->z, TotalAtomNumber*sizeof(float));
           cudaMalloc(&NucData_gpu->Nuc, TotalAtomNumber*sizeof(float));
           cudaMalloc(&NucData_gpu->K, sizeof(unsigned long int));
           cudaMalloc(&NucData_gpu->TotalAtomNumber, sizeof(unsigned long int));
           cudaMalloc(&NucData_gpu->RotMat, 9*sizeof(float));
           cudaMalloc(&NucData_gpu->NumberOfElements, *NucData->K*sizeof(unsigned int));
           cudaMalloc(&NucData_gpu->N_cum, *NucData->K*sizeof(unsigned long int));
           cudaMalloc(&NucData_gpu->N_act, *NucData->K*sizeof(unsigned long int));
           cudaMalloc(&NucData_gpu->N_avg, sizeof(unsigned long int));
     
           cudaMemcpy(NucData_gpu->N_avg, NucData->N_avg, sizeof(unsigned long int), cudaMemcpyHostToDevice);
           cudaMemcpy(NucData_gpu->K, NucData->K, sizeof(unsigned long int), cudaMemcpyHostToDevice);
           cudaMemcpy(NucData_gpu->TotalAtomNumber, NucData->TotalAtomNumber, sizeof(unsigned long int), cudaMemcpyHostToDevice);
           cudaMemcpy(NucData_gpu->RotMat, NucData->RotMat, 9*sizeof(float), cudaMemcpyHostToDevice);
           cudaMemcpy(NucData_gpu->NumberOfElements, NucData->NumberOfElements, *NucData->K*sizeof(unsigned int), cudaMemcpyHostToDevice);
           cudaMemcpy(NucData_gpu->N_cum, NucData->N_cum, *NucData->K*sizeof(unsigned long int), cudaMemcpyHostToDevice);
           cudaMemcpy(NucData_gpu->N_act, NucData->N_act, *NucData->K*sizeof(unsigned long int), cudaMemcpyHostToDevice);
     
       LogSystem::write("");
               // copy data from Host to Device
       LogSystem::write("copy NucData from RAM to GPU...");
           cudaMemcpy(NucData_gpu->x, NucData->x, TotalAtomNumber*sizeof(float), cudaMemcpyHostToDevice);
       LogSystem::write("  x done...");
           cudaMemcpy(NucData_gpu->y, NucData->y, TotalAtomNumber*sizeof(float), cudaMemcpyHostToDevice);
       LogSystem::write("  y done...");
           cudaMemcpy(NucData_gpu->z, NucData->z, TotalAtomNumber*sizeof(float), cudaMemcpyHostToDevice);
       LogSystem::write("  z done...");
           cudaMemcpy(NucData_gpu->Nuc, NucData->Nuc, TotalAtomNumber*sizeof(float), cudaMemcpyHostToDevice);
       LogSystem::write("  n done...");
       LogSystem::write("");
       LogSystem::write("data loaded...");
       LogSystem::write("");
    
          cudaMalloc(&NucData_gpu, sizeof(NuclearData));
          cudaMemcpy(NucData_gpu, NucData, sizeof(NuclearData), cudaMemcpyHostToDevice);
     
   }
   
   void read_NuclearData(NuclearData* NucData, \
                         NucDataProperties* NucDataProp, \
                         InputFileData* InputData, \
                         int NucData_File_Index){
   
       LogSystem::write("");
       LogSystem::write("load NucData...");
   
         unsigned long int K = *NucData->K;
   
         for(int i = 0; i < 9; i++){
            NucData->RotMat[i] = InputData->RotMat[i];
         }
   
         string filename;
         unsigned long int n = 0;
         float x_buf, y_buf, z_buf, Nuc_buf;
         ifstream fin;
         float x_mean = 0.0;
         float y_mean = 0.0;
         float z_mean = 0.0;
   
         unsigned long int N_act = 0;
         unsigned long int N_cum = 0;
   
   
         for(unsigned long int k = 1; k <= K; k++){
   
             N_act = NucData->NumberOfElements[k-1];
             NucData->N_act[k-1] = N_act;
   
             filename = NucDataProp->GlobalFolderPath + "/" + NucDataProp->SubFolderNames_Nom + "_" + to_string(k) \
                      + "/" + NucDataProp->SubFolder_FileNames_Nom + "_" + to_string(NucData_File_Index)\
                      + "." + NucDataProp->SubFolder_FileNames_Type;
   
       LogSystem::write(filename);
   
             fin.open(filename);
             n = 0;
             x_mean = 0.0;
             y_mean = 0.0;
             z_mean = 0.0;
   
            // read in the data and calculate the spatial average position
             while(fin >> x_buf >> y_buf >> z_buf >> Nuc_buf){
   
                     NucData->x[n + N_cum] = x_buf * InputData->XYZ_Unit_Factor;
                     NucData->y[n + N_cum] = y_buf * InputData->XYZ_Unit_Factor;
                     NucData->z[n + N_cum] = z_buf * InputData->XYZ_Unit_Factor;
                     NucData->Nuc[n + N_cum] = Nuc_buf;
      
                     x_mean += NucData->x[n + N_cum]/N_act;
                     y_mean += NucData->y[n + N_cum]/N_act;
                     z_mean += NucData->z[n + N_cum]/N_act;
   
                     n += 1;
             }
             fin.close();
   
             // centering of the xyz-data set
             for(int l = 0; l < N_act; l++){
                 NucData->x[l + N_cum] = NucData->x[l + N_cum] - x_mean;
                 NucData->y[l + N_cum] = NucData->y[l + N_cum] - y_mean;
                 NucData->z[l + N_cum] = NucData->z[l + N_cum] - z_mean;
             }
   
             N_cum += N_act;
             if(k<K){
                 NucData->N_cum[k] = N_cum;
             }
       LogSystem::write("N_act: " + std::to_string(N_act) + ", " + "N_cum: " + std::to_string(N_cum));
          }
   
         *NucData->N_avg = (unsigned long int) (((float)N_cum)/((float) K));
       LogSystem::write("N_avg: " + std::to_string(*NucData->N_avg));
       LogSystem::write("read (x,y,z,n) data finished...");
   
   }
   
   
   
   
   void init_NuclearData(NuclearData* NucData, \
                         NuclearData* NucData_gpu, \
                         NucDataProperties* NucDataProp, \
                         InputFileData* InputData, \
                         int NucData_File_Index){
   
        allocate_NuclearDataRAM(NucData, NucDataProp, NucData_File_Index);
        read_NuclearData(NucData, NucDataProp, InputData, NucData_File_Index);
        allocate_NuclearDataGPU(NucData, NucData_gpu, NucData_File_Index);
   }
   
   void free_NuclearData(NuclearData *NucData, \
                         NuclearData *NucData_gpu){
   
       LogSystem::write("free NucData...");
   
   //  cudaError_t err;
   
         free(NucData->x);
         free(NucData->y);
         free(NucData->z);
         free(NucData->Nuc);
         free(NucData->K);
         free(NucData->TotalAtomNumber);
         free(NucData->RotMat);
         free(NucData->NumberOfElements);
         free(NucData->N_cum);
         free(NucData->N_act);
         free(NucData->N_avg);
   
         cudaDeviceSynchronize();
         cudaFree(NucData_gpu->x);
         cudaDeviceSynchronize();
         cudaFree(NucData_gpu->y);
         cudaDeviceSynchronize();
         cudaFree(NucData_gpu->z);
         cudaDeviceSynchronize();
         cudaFree(NucData_gpu->Nuc);
         cudaDeviceSynchronize();
         cudaFree(NucData_gpu->K);
         cudaDeviceSynchronize();
         cudaFree(NucData_gpu->TotalAtomNumber);
         cudaDeviceSynchronize();
         cudaFree(NucData_gpu->RotMat);
         cudaDeviceSynchronize();
         cudaFree(NucData_gpu->NumberOfElements);
         cudaDeviceSynchronize();
         cudaFree(NucData_gpu->N_cum);
         cudaDeviceSynchronize();
         cudaFree(NucData_gpu->N_act);
         cudaDeviceSynchronize();
         cudaFree(NucData_gpu->N_avg);
         cudaDeviceSynchronize();
   
       LogSystem::write("free NucData finished.");
   
   }
   
   void disp_NuclearData(NuclearData *NucData){
        for(int k=0; k < *NucData->TotalAtomNumber; k++){
            cout << NucData->x[k] << " "\
                 << NucData->y[k] << " "\
                 << NucData->z[k] << " "\
                 << NucData->Nuc[k] << " "\
                 << "\n";
        }
   }
   
   
