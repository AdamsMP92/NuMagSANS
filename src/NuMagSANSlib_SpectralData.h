// File         : NuMagSANSlib_SANSData.h
// Author       : Dr. Michael Philipp ADAMS
// Company      : University of Luxembourg
// Department   : Department of Physics and Materials Sciences
// Group        : NanoMagnetism Group
// Group Leader : Prof. Andreas Michels
// Version      : 28 October 2025
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

struct SpectralData {
   unsigned int *Nq;        // number of q-values
   unsigned int *Ntheta;    // number of theta-values
   unsigned int *k_max;     // maximum number of modes
   float *dtheta;           // angular step-size

   float *q;

   // the arrays are single column and contain the angular sine and cosine spectra
   // the total array length is 2 * Nq * k_max
   float *I_Nuc_unpolarized;		  // nuclear SANS cross section
  	float *I_Mag_unpolarized;		  // unpolarized magnetic SANS cross section
  	float *I_Mag_polarized;			  // polarized magnetic SANS cross section
  	float *I_NucMag;					      // nuclear-magnetic interference SANS cross section
  	float *I_Mag_spin_flip;			  // spin-flip magnetic SANS cross section
  	float *I_Mag_chiral;				    // chiral magnetic SANS cross section
  	float *I_Mag_spin_flip_pm;		  // pm-spin-flip magnetic SANS cross section
  	float *I_Mag_spin_flip_mp;		  // mp-spin-flip magnetic SANS cross section
  	float *I_Mag_non_spin_flip_pp;	// pp-non-spin-flip magnetic SANS cross section
  	float *I_Mag_non_spin_flip_mm;	// mm-non-spin-flip magnetic SANS cross section
  	float *I_Mag_sanspol_p;			  // p-sanspol magnetic SANS cross section
  	float *I_Mag_sanspol_m;			  // m-sanspol magnetic SANS cross section

};

 
void allocate_SpectralData_RAM(InputFileData* InputData, \
                               ScatteringData* SANSData,\ 
                               SpectralData* SpecData){

    SpecData->Nq = (unsigned int*) malloc(sizeof(unsigned int));
    SpecData->Ntheta = (unsigned int*) malloc(sizeof(unsigned int));
    SpecData->k_max = (unsigned int*) malloc(sizeof(unsigned int));
    SpecData->dtheta = (float*) malloc(sizeof(float));
 
    *SpecData->Nq = InputData->Nq;
    *SpecData->Ntheta = InputData->Ntheta;
    *SpecData->k_max = InputData->k_max;
    *SpecData->dtheta = *SANSData->dtheta;

    // allocate the q-space vector
    SpecData->q      = (float*) malloc((*(SpecData->Nq))*sizeof(float));
    memset(SpecData->q, 0, (*SpecData->Nq) * sizeof(float));
    memcpy(SpecData->q, SANSData->q_1D, (*SpecData->Nq) * sizeof(float));

    // allocate spectral arrays
    size_t len = 2 * (*(SpecData->Nq)) * ((*(SpecData->k_max)) + 1) * sizeof(float);
    SpecData->I_Nuc_unpolarized      = (float*) malloc(len);
    SpecData->I_Mag_unpolarized      = (float*) malloc(len);
    SpecData->I_Mag_polarized        = (float*) malloc(len);
    SpecData->I_NucMag               = (float*) malloc(len);
    SpecData->I_Mag_spin_flip        = (float*) malloc(len);
    SpecData->I_Mag_chiral           = (float*) malloc(len);
    SpecData->I_Mag_spin_flip_pm     = (float*) malloc(len);
    SpecData->I_Mag_spin_flip_mp     = (float*) malloc(len);
    SpecData->I_Mag_non_spin_flip_pp = (float*) malloc(len);
    SpecData->I_Mag_non_spin_flip_mm = (float*) malloc(len);
    SpecData->I_Mag_sanspol_p        = (float*) malloc(len);
    SpecData->I_Mag_sanspol_m        = (float*) malloc(len);

    // --- Initialize to zero ---
    memset(SpecData->I_Nuc_unpolarized, 0, len);
    memset(SpecData->I_Mag_unpolarized, 0, len);
    memset(SpecData->I_Mag_polarized, 0, len);
    memset(SpecData->I_NucMag, 0, len);
    memset(SpecData->I_Mag_spin_flip, 0, len);
    memset(SpecData->I_Mag_chiral, 0, len);
    memset(SpecData->I_Mag_spin_flip_pm, 0, len);
    memset(SpecData->I_Mag_spin_flip_mp, 0, len);
    memset(SpecData->I_Mag_non_spin_flip_pp, 0, len);
    memset(SpecData->I_Mag_non_spin_flip_mm, 0, len);
    memset(SpecData->I_Mag_sanspol_p, 0, len);
    memset(SpecData->I_Mag_sanspol_m, 0, len);
 
}

void allocate_SpectralData_GPU(SpectralData* SpecData, \
                               SpectralData* SpecData_gpu){
 
    // Allocate scalar members on GPU
    cudaMalloc(&(SpecData_gpu->Nq),     sizeof(unsigned int));
    cudaMalloc(&(SpecData_gpu->Ntheta), sizeof(unsigned int));
    cudaMalloc(&(SpecData_gpu->k_max),  sizeof(unsigned int));
    cudaMalloc(&(SpecData_gpu->dtheta), sizeof(float));

    // Copy scalar values from host SpecData → device SpecData_gpu
    cudaMemcpy(SpecData_gpu->Nq,     SpecData->Nq,     sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMemcpy(SpecData_gpu->Ntheta, SpecData->Ntheta, sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMemcpy(SpecData_gpu->k_max,  SpecData->k_max,  sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMemcpy(SpecData_gpu->dtheta, SpecData->dtheta, sizeof(float),        cudaMemcpyHostToDevice);

    unsigned int Nq    = *SpecData->Nq;
    unsigned int k_max = *SpecData->k_max;
    size_t len_q  = Nq * sizeof(float);
    size_t len_I  = 2 * Nq * (k_max + 1) * sizeof(float);

    // Allocate arrays on GPU
    cudaMalloc((void**)&(SpecData_gpu->q), len_q);
    cudaMemcpy(SpecData_gpu->q, SpecData->q, len_q, cudaMemcpyHostToDevice);

    cudaMalloc(&(SpecData_gpu->I_Nuc_unpolarized), len_I);
    cudaMalloc(&(SpecData_gpu->I_Mag_unpolarized), len_I);
    cudaMalloc(&(SpecData_gpu->I_Mag_polarized), len_I);
    cudaMalloc(&(SpecData_gpu->I_NucMag), len_I);
    cudaMalloc(&(SpecData_gpu->I_Mag_spin_flip), len_I);
    cudaMalloc(&(SpecData_gpu->I_Mag_chiral), len_I);
    cudaMalloc(&(SpecData_gpu->I_Mag_spin_flip_pm), len_I);
    cudaMalloc(&(SpecData_gpu->I_Mag_spin_flip_mp), len_I);
    cudaMalloc(&(SpecData_gpu->I_Mag_non_spin_flip_pp), len_I);
    cudaMalloc(&(SpecData_gpu->I_Mag_non_spin_flip_mm), len_I);
    cudaMalloc(&(SpecData_gpu->I_Mag_sanspol_p), len_I);
    cudaMalloc(&(SpecData_gpu->I_Mag_sanspol_m), len_I);

    // Initialize GPU arrays with zeros
    cudaMemset(SpecData_gpu->I_Nuc_unpolarized,      0, len_I);
    cudaMemset(SpecData_gpu->I_Mag_unpolarized,      0, len_I);
    cudaMemset(SpecData_gpu->I_Mag_polarized,        0, len_I);
    cudaMemset(SpecData_gpu->I_NucMag,               0, len_I);
    cudaMemset(SpecData_gpu->I_Mag_spin_flip,        0, len_I);
    cudaMemset(SpecData_gpu->I_Mag_chiral,           0, len_I);
    cudaMemset(SpecData_gpu->I_Mag_spin_flip_pm,     0, len_I);
    cudaMemset(SpecData_gpu->I_Mag_spin_flip_mp,     0, len_I);
    cudaMemset(SpecData_gpu->I_Mag_non_spin_flip_pp, 0, len_I);
    cudaMemset(SpecData_gpu->I_Mag_non_spin_flip_mm, 0, len_I);
    cudaMemset(SpecData_gpu->I_Mag_sanspol_p,        0, len_I);
    cudaMemset(SpecData_gpu->I_Mag_sanspol_m,        0, len_I);

    cudaMalloc(&SpecData_gpu, sizeof(SpectralData));
    cudaMemcpy(SpecData_gpu, SpecData, sizeof(ScatteringData), cudaMemcpyHostToDevice);
 
}

void init_SpectralData(InputFileData *InputData, \
                       ScatteringData* SANSData, \
                       SpectralData* SpecData, \
                       SpectralData* SpecData_gpu){

     allocate_SpectralData_RAM(InputData, SANSData, SpecData);
     allocate_SpectralData_GPU(SpecData, SpecData_gpu);

}


void copyGPU2RAM_SpectralData(SpectralData *S, \
							                       SpectralData *S_gpu){

    cudaDeviceSynchronize();

    unsigned int Nq    = *S->Nq;
    unsigned int k_max = *S->k_max;

    size_t len_I = 2 * Nq * (k_max + 1) * sizeof(float);

    // Copy spectral arrays 
    cudaMemcpy(S->I_Nuc_unpolarized,      S_gpu->I_Nuc_unpolarized,      len_I, cudaMemcpyDeviceToHost);
    cudaMemcpy(S->I_Mag_unpolarized,      S_gpu->I_Mag_unpolarized,      len_I, cudaMemcpyDeviceToHost);
    cudaMemcpy(S->I_Mag_polarized,        S_gpu->I_Mag_polarized,        len_I, cudaMemcpyDeviceToHost);
    cudaMemcpy(S->I_NucMag,               S_gpu->I_NucMag,               len_I, cudaMemcpyDeviceToHost);
    cudaMemcpy(S->I_Mag_spin_flip,        S_gpu->I_Mag_spin_flip,        len_I, cudaMemcpyDeviceToHost);
    cudaMemcpy(S->I_Mag_chiral,           S_gpu->I_Mag_chiral,           len_I, cudaMemcpyDeviceToHost);
    cudaMemcpy(S->I_Mag_spin_flip_pm,     S_gpu->I_Mag_spin_flip_pm,     len_I, cudaMemcpyDeviceToHost);
    cudaMemcpy(S->I_Mag_spin_flip_mp,     S_gpu->I_Mag_spin_flip_mp,     len_I, cudaMemcpyDeviceToHost);
    cudaMemcpy(S->I_Mag_non_spin_flip_pp, S_gpu->I_Mag_non_spin_flip_pp, len_I, cudaMemcpyDeviceToHost);
    cudaMemcpy(S->I_Mag_non_spin_flip_mm, S_gpu->I_Mag_non_spin_flip_mm, len_I, cudaMemcpyDeviceToHost);
    cudaMemcpy(S->I_Mag_sanspol_p,        S_gpu->I_Mag_sanspol_p,        len_I, cudaMemcpyDeviceToHost);
    cudaMemcpy(S->I_Mag_sanspol_m,        S_gpu->I_Mag_sanspol_m,        len_I, cudaMemcpyDeviceToHost);

    // === 4. Synchronize to ensure data is ready ===
    cudaDeviceSynchronize();

}


void write2CSVtable_SpectralData(InputFileData *InputData, \
					     	     SpectralData *SpecData, \
					     	     int MagData_File_Index){


}



void free_SpectralData(SpectralData* S, \
                       SpectralData* S_gpu)
{

    free(S->Nq);
    free(S->Ntheta);
    free(S->k_max);
    free(S->dtheta);
    free(S->q);
    free(S->I_Nuc_unpolarized);
    free(S->I_Mag_unpolarized);
    free(S->I_Mag_polarized);
    free(S->I_NucMag);
    free(S->I_Mag_spin_flip);
    free(S->I_Mag_chiral);
    free(S->I_Mag_spin_flip_pm);
    free(S->I_Mag_spin_flip_mp);
    free(S->I_Mag_non_spin_flip_pp);
    free(S->I_Mag_non_spin_flip_mm);
    free(S->I_Mag_sanspol_p);
    free(S->I_Mag_sanspol_m);

    cudaDeviceSynchronize();
    cudaFree(S_gpu->Nq); 
    cudaDeviceSynchronize();
    cudaFree(S_gpu->Ntheta); 
    cudaDeviceSynchronize();
    cudaFree(S_gpu->k_max);
    cudaDeviceSynchronize();
    cudaFree(S_gpu->dtheta); 
    cudaDeviceSynchronize();
    cudaFree(S_gpu->q);
    cudaDeviceSynchronize();
    cudaFree(S_gpu->I_Nuc_unpolarized); 
    cudaDeviceSynchronize();
    cudaFree(S_gpu->I_Mag_unpolarized);
    cudaDeviceSynchronize();
    cudaFree(S_gpu->I_Mag_polarized); 
    cudaDeviceSynchronize();
    cudaFree(S_gpu->I_NucMag);
    cudaDeviceSynchronize();
    cudaFree(S_gpu->I_Mag_spin_flip); 
    cudaDeviceSynchronize();
    cudaFree(S_gpu->I_Mag_chiral);
    cudaDeviceSynchronize();
    cudaFree(S_gpu->I_Mag_spin_flip_pm); 
    cudaDeviceSynchronize();
    cudaFree(S_gpu->I_Mag_spin_flip_mp);
    cudaDeviceSynchronize();
    cudaFree(S_gpu->I_Mag_non_spin_flip_pp); 
    cudaDeviceSynchronize();
    cudaFree(S_gpu->I_Mag_non_spin_flip_mm);
    cudaDeviceSynchronize();
    cudaFree(S_gpu->I_Mag_sanspol_p); 
    cudaDeviceSynchronize();
    cudaFree(S_gpu->I_Mag_sanspol_m);
    cudaDeviceSynchronize();

}



