// File         : NuMagSANSlib_MagData.h
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


struct MagnetizationData{

   float *RotMat;	// rotation matrix R = [R0, R3, R6; R1, R4, R7; R2, R5, R8]
   float*x;		// position data x
   float*y;		// position data y
   float*z;		// position data z
   float*mx;	// magnetization data mx
   float*my;	// magnetization data my
   float*mz;	// magnetization data mz

   unsigned long int* K; 	// number of objects (e.g., orientations)
   unsigned long int* TotalAtomNumber;	// total atom number

   unsigned int* NumberOfElements;
   unsigned int* NumberOfNonZeroMoments;
   unsigned long int* N_act;    // number of atoms in k-th object
   unsigned long int* N_cum;    // cummulated number of atoms from object to object
   unsigned long int* N_avg;    // average number of atoms per object
};





void allocate_MagnetizationDataRAM(MagnetizationData* MagData, \
						           MagDataProperties* MagDataProp, \
						           InputFileData* InputData, \
                                   int MagData_File_Index){

	 bool ExcludeZeroMoments_flag = InputData->ExcludeZeroMoments_flag;
	 unsigned long int K = MagDataProp->Number_Of_SubFolders;
     unsigned long int TotalAtomNumber = 0;
	 if(ExcludeZeroMoments_flag){
     	TotalAtomNumber = MagDataProp->TotalNZMAtomNumber[MagData_File_Index];
     } else{
		TotalAtomNumber = MagDataProp->TotalAtomNumber[MagData_File_Index];
     }

     MagData->x = (float*) malloc(TotalAtomNumber*sizeof(float));
     MagData->y = (float*) malloc(TotalAtomNumber*sizeof(float));
     MagData->z = (float*) malloc(TotalAtomNumber*sizeof(float));
     MagData->mx = (float*) malloc(TotalAtomNumber*sizeof(float));
     MagData->my = (float*) malloc(TotalAtomNumber*sizeof(float));
     MagData->mz = (float*) malloc(TotalAtomNumber*sizeof(float));
     MagData->K = (unsigned long int*) malloc(sizeof(unsigned long int));
     MagData->TotalAtomNumber = (unsigned long int*) malloc(sizeof(unsigned long int));
     MagData->RotMat = (float*) malloc(9*sizeof(float));
     MagData->NumberOfElements = (unsigned int*) malloc(K*sizeof(unsigned int));
     MagData->NumberOfNonZeroMoments = (unsigned int*) malloc(K*sizeof(unsigned int));
     MagData->N_cum = (unsigned long int*) malloc(K*sizeof(unsigned long int));
     MagData->N_act = (unsigned long int*) malloc(K*sizeof(unsigned long int));
     MagData->N_avg = (unsigned long int*) malloc(sizeof(unsigned long int));

	 if (!MagData->x || !MagData->y || !MagData->z || 
	     !MagData->mx || !MagData->my || !MagData->mz || 
	     !MagData->N_avg || !MagData->K || !MagData->TotalAtomNumber || !MagData->RotMat) {
	         perror("Memory allocation failed");
	         exit(EXIT_FAILURE);
	 }

	 *MagData->N_avg = 0;
	 *MagData->K = K;
	 *MagData->TotalAtomNumber = TotalAtomNumber;

	 for(unsigned long int i = 0; i < TotalAtomNumber; i++){
		MagData->x[i] = 0.0;
		MagData->y[i] = 0.0;
		MagData->z[i] = 0.0;
		MagData->mx[i] = 0.0;
		MagData->my[i] = 0.0;
		MagData->mz[i] = 0.0;
	 }	

	 for(unsigned int i = 0; i < K; i++){
         MagData->NumberOfElements[i] = MagDataProp->NumberOfElements[i][MagData_File_Index];
         MagData->NumberOfNonZeroMoments[i] = MagDataProp->NumberOfNonZeroMoments[i][MagData_File_Index];
         MagData->N_cum[i] = 0;
         MagData->N_act[i] = 0;
    }

}


void allocate_MagnetizationDataGPU(MagnetizationData* MagData, \
 						           MagnetizationData* MagData_gpu, \
                                   int MagData_File_Index){

     unsigned long int TotalAtomNumber = *MagData->TotalAtomNumber;
	
	 cudaMalloc(&MagData_gpu->x, TotalAtomNumber*sizeof(float));
     cudaMalloc(&MagData_gpu->y, TotalAtomNumber*sizeof(float));
     cudaMalloc(&MagData_gpu->z, TotalAtomNumber*sizeof(float));
     cudaMalloc(&MagData_gpu->mx, TotalAtomNumber*sizeof(float));
     cudaMalloc(&MagData_gpu->my, TotalAtomNumber*sizeof(float));
     cudaMalloc(&MagData_gpu->mz, TotalAtomNumber*sizeof(float));
     cudaMalloc(&MagData_gpu->K, sizeof(unsigned long int));
	 cudaMalloc(&MagData_gpu->TotalAtomNumber, sizeof(unsigned long int));
	 cudaMalloc(&MagData_gpu->RotMat, 9*sizeof(float));
     cudaMalloc(&MagData_gpu->NumberOfElements, *MagData->K*sizeof(unsigned int));
     cudaMalloc(&MagData_gpu->NumberOfNonZeroMoments, *MagData->K*sizeof(unsigned int));
     cudaMalloc(&MagData_gpu->N_cum, *MagData->K*sizeof(unsigned long int));
     cudaMalloc(&MagData_gpu->N_act, *MagData->K*sizeof(unsigned long int));
      cudaMalloc(&MagData_gpu->N_avg, sizeof(unsigned long int));

	 cudaMemcpy(MagData_gpu->N_avg, MagData->N_avg, sizeof(unsigned long int), cudaMemcpyHostToDevice);
	 cudaMemcpy(MagData_gpu->K, MagData->K, sizeof(unsigned long int), cudaMemcpyHostToDevice);
	 cudaMemcpy(MagData_gpu->TotalAtomNumber, MagData->TotalAtomNumber, sizeof(unsigned long int), cudaMemcpyHostToDevice);
	 cudaMemcpy(MagData_gpu->RotMat, MagData->RotMat, 9*sizeof(float), cudaMemcpyHostToDevice);
     cudaMemcpy(MagData_gpu->NumberOfElements, MagData->NumberOfElements, *MagData->K*sizeof(unsigned int), cudaMemcpyHostToDevice);
     cudaMemcpy(MagData_gpu->NumberOfElements, MagData->NumberOfElements, *MagData->K*sizeof(unsigned int), cudaMemcpyHostToDevice);
     cudaMemcpy(MagData_gpu->N_cum, MagData->N_cum, *MagData->K*sizeof(unsigned long int), cudaMemcpyHostToDevice);
     cudaMemcpy(MagData_gpu->N_act, MagData->N_act, *MagData->K*sizeof(unsigned long int), cudaMemcpyHostToDevice);
		
     cout << " \n";
         // copy data from Host to Device
     cout << "Copy Data from RAM to GPU...\n";
     cudaMemcpy(MagData_gpu->x, MagData->x, TotalAtomNumber*sizeof(float), cudaMemcpyHostToDevice);
     cout << "   x done...\n";
     cudaMemcpy(MagData_gpu->y, MagData->y, TotalAtomNumber*sizeof(float), cudaMemcpyHostToDevice);
     cout << "   y done...\n";
     cudaMemcpy(MagData_gpu->z, MagData->z, TotalAtomNumber*sizeof(float), cudaMemcpyHostToDevice);
     cout << "   z done...\n";
     cudaMemcpy(MagData_gpu->mx, MagData->mx, TotalAtomNumber*sizeof(float), cudaMemcpyHostToDevice);
     cout << "   mx done...\n";
     cudaMemcpy(MagData_gpu->my, MagData->my, TotalAtomNumber*sizeof(float), cudaMemcpyHostToDevice);
     cout << "   my done...\n";
     cudaMemcpy(MagData_gpu->mz, MagData->mz, TotalAtomNumber*sizeof(float), cudaMemcpyHostToDevice);
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
     cout << "Load Magnetization Data...\n";

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

     unsigned long int N_cum = 0;   // cummulative counter
     unsigned long int N_act = 0;   // actual number of atoms per file


     for(unsigned long int k = 1; k <= K; k++){

         if(ExcludeZeroMoments_flag){
             N_act = MagData->NumberOfNonZeroMoments[k-1];
         }else{
             N_act = MagData->NumberOfElements[k-1];
         }
         MagData->N_act[k-1] = N_act;

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
                	MagData->x[n + N_cum] = x_buf * InputData->XYZ_Unit_Factor;
                	MagData->y[n + N_cum] = y_buf * InputData->XYZ_Unit_Factor;
                	MagData->z[n + N_cum] = z_buf * InputData->XYZ_Unit_Factor;
                	MagData->mx[n + N_cum] = mx_buf;
                	MagData->my[n + N_cum] = my_buf;
                	MagData->mz[n + N_cum] = mz_buf;

                	x_mean += MagData->x[n + N_cum]/N_act;
                	y_mean += MagData->y[n + N_cum]/N_act;
                	z_mean += MagData->z[n + N_cum]/N_act;

                	n += 1;
             	}
         	}
         } else{
			while(fin >> x_buf >> y_buf >> z_buf >> mx_buf >> my_buf >> mz_buf){
			    
			        MagData->x[n + N_cum] = x_buf * InputData->XYZ_Unit_Factor;
			        MagData->y[n + N_cum] = y_buf * InputData->XYZ_Unit_Factor;
			        MagData->z[n + N_cum] = z_buf * InputData->XYZ_Unit_Factor;
			        MagData->mx[n + N_cum] = mx_buf;
			        MagData->my[n + N_cum] = my_buf;
			        MagData->mz[n + N_cum] = mz_buf;
			
			        x_mean += MagData->x[n + N_cum]/N_act;
			        y_mean += MagData->y[n + N_cum]/N_act;
			        z_mean += MagData->z[n + N_cum]/N_act;
			
                	n += 1;
			}        	
         }
         
         fin.close();

         // centering of the xyz-data set
         for(int l = 0; l < N_act; l++){
             MagData->x[l + N_cum] = MagData->x[l + N_cum] - x_mean;
             MagData->y[l + N_cum] = MagData->y[l + N_cum] - y_mean;
             MagData->z[l + N_cum] = MagData->z[l + N_cum] - z_mean;
         }
 
        // update of the cummulative counter
        N_cum += N_act;
        if(k<K){
            MagData->N_cum[k] = N_cum;
        }
        cout << "N_act: " << N_act << ",  " << "N_cum: " << N_cum << "\n";

      }

      *MagData->N_avg = (unsigned long int) (((float)N_cum)/((float) K));
      cout << "N_avg: " << *MagData->N_avg << "\n";

     cout << "(x, y, z, mx, my, mz) - data load done...\n";
    
 
 }

void init_MagnetizationData(MagnetizationData* MagData, \
                  MagnetizationData* MagData_gpu, \
                  MagDataProperties* MagDataProp, \
                  InputFileData* InputData, \
                  int MagData_File_Index){

	allocate_MagnetizationDataRAM(MagData, MagDataProp, InputData, MagData_File_Index);
	read_MagnetizationData(MagData, MagDataProp, InputData, MagData_File_Index);
	allocate_MagnetizationDataGPU(MagData, MagData_gpu, MagData_File_Index);
   
}

void free_MagnetizationData(MagnetizationData *MagData, \
                  MagnetizationData *MagData_gpu){
 

     cout << "Free Magnetization Data..." << "\n";

     free(MagData->x);
     free(MagData->y);
     free(MagData->z);
     free(MagData->mx);
     free(MagData->my);
     free(MagData->mz);
     free(MagData->K);
     free(MagData->TotalAtomNumber);
     free(MagData->RotMat);
     free(MagData->NumberOfElements);
     free(MagData->NumberOfNonZeroMoments);
     free(MagData->N_cum);
     free(MagData->N_act);
     free(MagData->N_avg);

     cudaFree(MagData_gpu->x);
     cudaFree(MagData_gpu->y);
     cudaFree(MagData_gpu->z);
     cudaFree(MagData_gpu->mx);
     cudaFree(MagData_gpu->my);
     cudaFree(MagData_gpu->mz);
	 cudaFree(MagData_gpu->K);
	 cudaFree(MagData_gpu->TotalAtomNumber);
	 cudaFree(MagData_gpu->RotMat);
     cudaFree(MagData_gpu->NumberOfElements);
     cudaFree(MagData_gpu->NumberOfNonZeroMoments);
     cudaFree(MagData_gpu->N_cum);
     cudaFree(MagData_gpu->N_act);
     cudaFree(MagData_gpu->N_avg);

	 cudaFree(MagData_gpu);

}

void disp_MagnetizationData(MagnetizationData *MagData){
	for(int k=0; k < *MagData->TotalAtomNumber; k++){
		cout << MagData->x[k] << " "\
		     << MagData->y[k] << " "\
		     << MagData->z[k] << " "\
		     << MagData->mx[k] << " "\
		     << MagData->my[k] << " "\
		     << MagData->mz[k] << "\n";
	}
}
