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
    size_t len = 2 * (*(SpecData->Nq)) * (*(SpecData->k_max)) * sizeof(float);
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

 
}

void init_SpectralData(InputFileData *InputData, \
                       ScatteringData* SANSData, \
                       SpectralData* SpecData, \
                       SpectralData* SpecData_gpu){

     allocate_SpectralData_RAM(InputData, SANSData, SpecData);
     allocate_SpectralData_GPU(SpecData, SpecData_gpu;

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
 
}



