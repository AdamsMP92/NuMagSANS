// File         : NuMagSANSlib_NucData.h
// Author       : Michael Philipp ADAMS, M.Sc.
// Company      : University of Luxembourg
// Department   : Department of Physics and Materials Sciences
// Group        : NanoMagnetism Group
// Group Leader : Prof. Andreas Michels
// Version      : 22 October 2024
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
    unsigned long int* N;    // number of nuclear points
    unsigned long int* K;    // number of orientations
    unsigned long int* N_k;  // product of N * K

};

void allocate_NuclearDataRAM(NuclearData* NucData, \
                             NucDataProperties* NucDataProp){
   
         unsigned long int K = NucDataProp->Number_Of_SubFolders;
         unsigned long int N = NucDataProp->NumberOfElements[0][0];
         unsigned long int N_k = N * K;
   
         NucData->x = (float*) malloc(N_k*sizeof(float));
         NucData->y = (float*) malloc(N_k*sizeof(float));
         NucData->z = (float*) malloc(N_k*sizeof(float));
         NucData->Nuc = (float*) malloc(N_k*sizeof(float));
         NucData->N = (unsigned long int*) malloc(sizeof(unsigned long int));
         NucData->K = (unsigned long int*) malloc(sizeof(unsigned long int));
         NucData->N_k = (unsigned long int*) malloc(sizeof(unsigned long int));
         NucData->RotMat = (float*) malloc(9*sizeof(float));

		if(!NucData->x || !NucData->y || !NucData->z || !NucData->Nuc ||
		   !NucData->N || !NucData->K || !NucData->N_k || !NucData->RotMat){
		   		perror("Memory alloation failed");
		   		exit(EXIT_FAILURE);
		}
   
         *NucData->N = N;
         *NucData->K = K;
         *NucData->N_k = N_k;

         for(unsigned long int i = 0; i < N_k; i++){
         	NucData->x[i] = 0.0;
         	NucData->y[i] = 0.0;
         	NucData->z[i] = 0.0;
         	NucData->Nuc[i] = 0.0;
         }
}


 void allocate_NuclearDataGPU(NuclearData* NucData, \
                              NuclearData* NucData_gpu){
  
        unsigned long int N_k = *NucData->N_k;
  
        cudaMalloc(&NucData_gpu->x, N_k*sizeof(float));
        cudaMalloc(&NucData_gpu->y, N_k*sizeof(float));
        cudaMalloc(&NucData_gpu->z, N_k*sizeof(float));
        cudaMalloc(&NucData_gpu->Nuc, N_k*sizeof(float));
        cudaMalloc(&NucData_gpu->N, sizeof(unsigned long int));
        cudaMalloc(&NucData_gpu->K, sizeof(unsigned long int));
        cudaMalloc(&NucData_gpu->N_k, sizeof(unsigned long int));
        cudaMalloc(&NucData_gpu->RotMat, 9*sizeof(float));
  
        cudaMemcpy(NucData_gpu->N, NucData->N, sizeof(unsigned long int), cudaMemcpyHostToDevice);
        cudaMemcpy(NucData_gpu->K, NucData->K, sizeof(unsigned long int), cudaMemcpyHostToDevice);
        cudaMemcpy(NucData_gpu->N_k, NucData->N_k, sizeof(unsigned long int), cudaMemcpyHostToDevice);
        cudaMemcpy(NucData_gpu->RotMat, NucData->RotMat, 9*sizeof(float), cudaMemcpyHostToDevice);
  
        cout << " \n";
            // copy data from Host to Device
        cout << "Copy Data from RAM to GPU...\n";
        cudaMemcpy(NucData_gpu->x, NucData->x, N_k*sizeof(float), cudaMemcpyHostToDevice);
        cout << "   x done...\n";
        cudaMemcpy(NucData_gpu->y, NucData->y, N_k*sizeof(float), cudaMemcpyHostToDevice);
        cout << "   y done...\n";
        cudaMemcpy(NucData_gpu->z, NucData->z, N_k*sizeof(float), cudaMemcpyHostToDevice);
        cout << "   z done...\n";
        cudaMemcpy(NucData_gpu->Nuc, NucData->Nuc, N_k*sizeof(float), cudaMemcpyHostToDevice);
        cout << "   Nuc done...\n";
        cout << " \n";
        cout << "Data loaded...\n";
        cout << " \n";
  
       cudaMalloc(&NucData_gpu, sizeof(NuclearData));
       cudaMemcpy(NucData_gpu, NucData, sizeof(NuclearData), cudaMemcpyHostToDevice);
  
}

void read_NuclearData(NuclearData* NucData, \
                      NucDataProperties* NucDataProp, \
                   	  InputFileData* InputData, \
                      int NucData_File_Index){

      cout << " \n";
      cout << "Load Data...\n";

      unsigned long int N = *NucData->N;
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

      for(unsigned long int k = 1; k <= K; k++){

          filename = NucDataProp->GlobalFolderPath + "/" + NucDataProp->SubFolderNames_Nom + "_" + to_string(k) \
                   + "/" + NucDataProp->SubFolder_FileNames_Nom + "_" + to_string(NucData_File_Index)\
                   + "." + NucDataProp->SubFolder_FileNames_Type;

          cout << filename << "\n";


          fin.open(filename);
          n = 0;
          x_mean = 0.0;
          y_mean = 0.0;
          z_mean = 0.0;

         // read in the data and calculate the spatial average position
          while(fin >> x_buf >> y_buf >> z_buf >> Nuc_buf){

                  NucData->x[n + (k-1)*N] = x_buf * InputData->XYZ_Unit_Factor;
                  NucData->y[n + (k-1)*N] = y_buf * InputData->XYZ_Unit_Factor;
                  NucData->z[n + (k-1)*N] = z_buf * InputData->XYZ_Unit_Factor;
                  NucData->Nuc[n + (k-1)*N] = Nuc_buf;
   
                  x_mean += NucData->x[n + (k-1)*N]/N;
                  y_mean += NucData->y[n + (k-1)*N]/N;
                  z_mean += NucData->z[n + (k-1)*N]/N;

                  n += 1;
          }
          fin.close();

          // centering of the xyz-data set
          for(int l = 0; l < N; l++){
              NucData->x[l + (k-1)*N] = NucData->x[l + (k-1)*N] - x_mean;
              NucData->y[l + (k-1)*N] = NucData->y[l + (k-1)*N] - y_mean;
              NucData->z[l + (k-1)*N] = NucData->z[l + (k-1)*N] - z_mean;
          }

       }

      cout << "(x, y, z, Nuc) - data load done...\n";
}




void init_NuclearData(NuclearData* NucData, \
                      NuclearData* NucData_gpu, \
                   	  NucDataProperties* NucDataProp, \
                  	  InputFileData* InputData, \
                   	  int NucData_File_Index){

     allocate_NuclearDataRAM(NucData, NucDataProp);
     read_NuclearData(NucData, NucDataProp, InputData, NucData_File_Index);
     allocate_NuclearDataGPU(NucData, NucData_gpu);
}

void free_NuclearData(NuclearData *NucData, \
                      NuclearData *NucData_gpu){

      free(NucData->x);
      free(NucData->y);
      free(NucData->z);
      free(NucData->Nuc);
      free(NucData->N);
      free(NucData->K);
      free(NucData->N_k);
      free(NucData->RotMat);

      cudaFree(NucData_gpu->x);
      cudaFree(NucData_gpu->y);
      cudaFree(NucData_gpu->z);
      cudaFree(NucData_gpu->Nuc);
      cudaFree(NucData_gpu->N);
      cudaFree(NucData_gpu->K);
      cudaFree(NucData_gpu->N_k);
      cudaFree(NucData_gpu->RotMat);

      cudaFree(NucData);

}

void disp_NuclearData(NuclearData *NucData){
     for(int k=0; k < *NucData->N_k; k++){
         cout << NucData->x[k] << " "\
              << NucData->y[k] << " "\
              << NucData->z[k] << " "\
              << NucData->Nuc[k] << " "\
              << "\n";
     }
}


