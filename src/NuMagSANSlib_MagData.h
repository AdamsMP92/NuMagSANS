// File         : NuMagSANSlib_MagData.h
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

struct MagnetizationData{

   float *RotMat;	// rotation matrix R = [R0, R3, R6; R1, R4, R7; R2, R5, R8]
   float*x;		// position data x
   float*y;		// position data y
   float*z;		// position data z
   float*mx;	// magnetization data mx
   float*my;	// magnetization data my
   float*mz;	// magnetization data mz
   unsigned long int* N;	// number of non-zero magnetic moments
   unsigned long int* K; 	// number of orientations
   unsigned long int* N_k;	// product of N * K
   
};

void allocate_MagnetizationDataRAM(MagnetizationData* MagData, \
						           MagDataProperties* MagDataProp, \
						           InputFileData* InputData){

	 bool ExcludeZeroMoments_flag = InputData->ExcludeZeroMoments_flag;
	 unsigned long int K = MagDataProp->Number_Of_SubFolders;
	 unsigned long int N;
	 if(ExcludeZeroMoments_flag){
     	N = MagDataProp->NumberOfNonZeroMoments[0][0];
     } else{
		N = MagDataProp->NumberOfElements[0][0];
     }
     unsigned long int N_k = N * K;
         
     MagData->x = (float*) malloc(N_k*sizeof(float));
     MagData->y = (float*) malloc(N_k*sizeof(float));
     MagData->z = (float*) malloc(N_k*sizeof(float));
     MagData->mx = (float*) malloc(N_k*sizeof(float));
     MagData->my = (float*) malloc(N_k*sizeof(float));
     MagData->mz = (float*) malloc(N_k*sizeof(float));
     MagData->N = (unsigned long int*) malloc(sizeof(unsigned long int));
     MagData->K = (unsigned long int*) malloc(sizeof(unsigned long int));
     MagData->N_k = (unsigned long int*) malloc(sizeof(unsigned long int));
     MagData->RotMat = (float*) malloc(9*sizeof(float));

	 if (!MagData->x || !MagData->y || !MagData->z || 
	     !MagData->mx || !MagData->my || !MagData->mz || 
	     !MagData->N || !MagData->K || !MagData->N_k || !MagData->RotMat) {
	         perror("Memory allocation failed");
	         exit(EXIT_FAILURE);
	 }

	 *MagData->N = N;
	 *MagData->K = K;
	 *MagData->N_k = N_k;

	 for(unsigned long int i = 0; i < N_k; i++){
		MagData->x[i] = 0.0;
		MagData->y[i] = 0.0;
		MagData->z[i] = 0.0;
		MagData->mx[i] = 0.0;
		MagData->my[i] = 0.0;
		MagData->mz[i] = 0.0;
	 }	
}


void allocate_MagnetizationDataGPU(MagnetizationData* MagData, \
 						           MagnetizationData* MagData_gpu){

     unsigned long int N_k = *MagData->N_k;
	
	 cudaMalloc(&MagData_gpu->x, N_k*sizeof(float));
     cudaMalloc(&MagData_gpu->y, N_k*sizeof(float));
     cudaMalloc(&MagData_gpu->z, N_k*sizeof(float));
     cudaMalloc(&MagData_gpu->mx, N_k*sizeof(float));
     cudaMalloc(&MagData_gpu->my, N_k*sizeof(float));
     cudaMalloc(&MagData_gpu->mz, N_k*sizeof(float));
     cudaMalloc(&MagData_gpu->N, sizeof(unsigned long int));
     cudaMalloc(&MagData_gpu->K, sizeof(unsigned long int));
	 cudaMalloc(&MagData_gpu->N_k, sizeof(unsigned long int));
	 cudaMalloc(&MagData_gpu->RotMat, 9*sizeof(float));

	 cudaMemcpy(MagData_gpu->N, MagData->N, sizeof(unsigned long int), cudaMemcpyHostToDevice);
	 cudaMemcpy(MagData_gpu->K, MagData->K, sizeof(unsigned long int), cudaMemcpyHostToDevice);
	 cudaMemcpy(MagData_gpu->N_k, MagData->N_k, sizeof(unsigned long int), cudaMemcpyHostToDevice);
	 cudaMemcpy(MagData_gpu->RotMat, MagData->RotMat, 9*sizeof(float), cudaMemcpyHostToDevice);
		
     cout << " \n";
         // copy data from Host to Device
     cout << "Copy Data from RAM to GPU...\n";
     cudaMemcpy(MagData_gpu->x, MagData->x, N_k*sizeof(float), cudaMemcpyHostToDevice);
     cout << "   x done...\n";
     cudaMemcpy(MagData_gpu->y, MagData->y, N_k*sizeof(float), cudaMemcpyHostToDevice);
     cout << "   y done...\n";
     cudaMemcpy(MagData_gpu->z, MagData->z, N_k*sizeof(float), cudaMemcpyHostToDevice);
     cout << "   z done...\n";
     cudaMemcpy(MagData_gpu->mx, MagData->mx, N_k*sizeof(float), cudaMemcpyHostToDevice);
     cout << "   mx done...\n";
     cudaMemcpy(MagData_gpu->my, MagData->my, N_k*sizeof(float), cudaMemcpyHostToDevice);
     cout << "   my done...\n";
     cudaMemcpy(MagData_gpu->mz, MagData->mz, N_k*sizeof(float), cudaMemcpyHostToDevice);
     cout << "   mz done...\n";
     cout << " \n";
     cout << "Data loaded...\n";
     cout << " \n";

	cudaMalloc(&MagData_gpu, sizeof(MagnetizationData)); 
	cudaMemcpy(MagData_gpu, MagData, sizeof(MagnetizationData), cudaMemcpyHostToDevice);
     
}


void read_MagnetizationData(MagnetizationData* MagData, \
				            MagDataProperties* MagDataProp, \
				            InputFileData* InputData, \
                            int MagData_File_Index){

	 cout << " \n";
     cout << "Load Data...\n";

     unsigned long int N = *MagData->N;
     unsigned long int K = *MagData->K;
     bool ExcludeZeroMoments_flag = InputData->ExcludeZeroMoments_flag;

	 for(int i = 0; i < 9; i++){
	 	MagData->RotMat[i] = InputData->RotMat[i];
	 }

     string filename;
     unsigned long int n = 0;
     float x_buf, y_buf, z_buf, mx_buf, my_buf, mz_buf;
     ifstream fin;
     float x_mean = 0.0;
     float y_mean = 0.0;
     float z_mean = 0.0;

     for(unsigned long int k = 1; k <= K; k++){

         filename = MagDataProp->GlobalFolderPath + "/" + MagDataProp->SubFolderNames_Nom + "_" + to_string(k) \
                  + "/" + MagDataProp->SubFolder_FileNames_Nom + "_" + to_string(MagData_File_Index)\
                  + "." + MagDataProp->SubFolder_FileNames_Type;
                  
         cout << filename << "\n";


         fin.open(filename);
         n = 0;
         x_mean = 0.0;
         y_mean = 0.0;
         z_mean = 0.0;

         // read in the data and calculate the spatial average position
         if(ExcludeZeroMoments_flag){
         	while(fin >> x_buf >> y_buf >> z_buf >> mx_buf >> my_buf >> mz_buf){
            	if(mx_buf != 0.0 || my_buf != 0.0 || mz_buf != 0.0){
                	MagData->x[n + (k-1)*N] = x_buf * InputData->XYZ_Unit_Factor;
                	MagData->y[n + (k-1)*N] = y_buf * InputData->XYZ_Unit_Factor;
                	MagData->z[n + (k-1)*N] = z_buf * InputData->XYZ_Unit_Factor;
                	MagData->mx[n + (k-1)*N] = mx_buf;
                	MagData->my[n + (k-1)*N] = my_buf;
                	MagData->mz[n + (k-1)*N] = mz_buf;

                	x_mean += MagData->x[n + (k-1)*N]/N;
                	y_mean += MagData->y[n + (k-1)*N]/N;
                	z_mean += MagData->z[n + (k-1)*N]/N;

                	n += 1;
             	}
         	}
         } else{
			while(fin >> x_buf >> y_buf >> z_buf >> mx_buf >> my_buf >> mz_buf){
			    
			        MagData->x[n + (k-1)*N] = x_buf * InputData->XYZ_Unit_Factor;
			        MagData->y[n + (k-1)*N] = y_buf * InputData->XYZ_Unit_Factor;
			        MagData->z[n + (k-1)*N] = z_buf * InputData->XYZ_Unit_Factor;
			        MagData->mx[n + (k-1)*N] = mx_buf;
			        MagData->my[n + (k-1)*N] = my_buf;
			        MagData->mz[n + (k-1)*N] = mz_buf;
			
			        x_mean += MagData->x[n + (k-1)*N]/N;
			        y_mean += MagData->y[n + (k-1)*N]/N;
			        z_mean += MagData->z[n + (k-1)*N]/N;
			
                	n += 1;
			}        	
         }
         
         fin.close();

         // centering of the xyz-data set
         for(int l = 0; l < N; l++){
             MagData->x[l + (k-1)*N] = MagData->x[l + (k-1)*N] - x_mean;
             MagData->y[l + (k-1)*N] = MagData->y[l + (k-1)*N] - y_mean;
             MagData->z[l + (k-1)*N] = MagData->z[l + (k-1)*N] - z_mean;
         }
 
      }

     cout << "(x, y, z, mx, my, mz) - data load done...\n";
    
 
 }

void init_MagnetizationData(MagnetizationData* MagData, \
                  MagnetizationData* MagData_gpu, \
                  MagDataProperties* MagDataProp, \
                  InputFileData* InputData, \
                  int MagData_File_Index){

	allocate_MagnetizationDataRAM(MagData, MagDataProp, InputData);
	read_MagnetizationData(MagData, MagDataProp, InputData, MagData_File_Index);
	allocate_MagnetizationDataGPU(MagData, MagData_gpu);
   
}

void free_MagnetizationData(MagnetizationData *MagData, \
                  MagnetizationData *MagData_gpu){
 
     free(MagData->x);
     free(MagData->y);
     free(MagData->z);
     free(MagData->mx);
     free(MagData->my);
     free(MagData->mz);
     free(MagData->N);
     free(MagData->K);
     free(MagData->N_k);
     free(MagData->RotMat);

     cudaFree(MagData_gpu->x);
     cudaFree(MagData_gpu->y);
     cudaFree(MagData_gpu->z);
     cudaFree(MagData_gpu->mx);
     cudaFree(MagData_gpu->my);
     cudaFree(MagData_gpu->mz);
	 cudaFree(MagData_gpu->N);
	 cudaFree(MagData_gpu->K);
	 cudaFree(MagData_gpu->N_k);
	 cudaFree(MagData_gpu->RotMat);

	 cudaFree(MagData);

}

void disp_MagnetizationData(MagnetizationData *MagData){
	for(int k=0; k < *MagData->N_k; k++){
		cout << MagData->x[k] << " "\
		     << MagData->y[k] << " "\
		     << MagData->z[k] << " "\
		     << MagData->mx[k] << " "\
		     << MagData->my[k] << " "\
		     << MagData->mz[k] << "\n";
	}
}
