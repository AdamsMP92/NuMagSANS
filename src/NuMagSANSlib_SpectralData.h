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
    unsigned int *Nq;        // Anzahl q-Werte
    unsigned int *Ntheta;    // Anzahl theta-Werte
    unsigned int *k_max;     // maximaler k-Wert
    float *dtheta; // Winkel-Schrittweite

    // the arrays are single column and contain the angular sine and cosine spectra
    // the total array length is 2 * Nq * k_max
    float *I_Nuc_2D_unpolarized;		  // nuclear SANS cross section
  	float *I_Mag_2D_unpolarized;		  // unpolarized magnetic SANS cross section
  	float *I_Mag_2D_polarized;			  // polarized magnetic SANS cross section
  	float *I_NucMag_2D;					      // nuclear-magnetic interference SANS cross section
  	float *I_Mag_2D_spin_flip;			  // spin-flip magnetic SANS cross section
  	float *I_Mag_2D_chiral;				    // chiral magnetic SANS cross section
  	float *I_Mag_2D_spin_flip_pm;		  // pm-spin-flip magnetic SANS cross section
  	float *I_Mag_2D_spin_flip_mp;		  // mp-spin-flip magnetic SANS cross section
  	float *I_Mag_2D_non_spin_flip_pp;	// pp-non-spin-flip magnetic SANS cross section
  	float *I_Mag_2D_non_spin_flip_mm;	// mm-non-spin-flip magnetic SANS cross section
  	float *I_Mag_2D_sanspol_p;			  // p-sanspol magnetic SANS cross section
  	float *I_Mag_2D_sanspol_m;			  // m-sanspol magnetic SANS cross section

};

 
void allocate_SpectralData_RAM(InputFileData* InputData, \
                               ScatteringData* SANSData,\ 
                               SpectralData* SpecData){

    SpectralData->Nq = (unsigned int*) malloc(sizeof(unsigned int));
    SpectralData->Ntheta = (unsigned int*) malloc(sizeof(unsigned int));
    SpectralData->k_max = (unsigned int*) malloc(sizeof(unsigned int));
    SpectralData->dtheta = (float*) malloc(sizeof(float));
 
    *SpectralData->Nq = InputData->Nq;
    *SpectralData->Ntheta = InputData->Ntheta;
    *SpectralData->k_max = InputData->k_max;
    *SpectralData->dtheta = *SANSData->dtheta;

    // --- Allocate spectral arrays ---
    size_t len = 2 * (*(SpecData->Nq)) * (*(SpecData->k_max)) * sizeof(float);

    SpecData->I_Nuc_2D_unpolarized      = (float*) malloc(len);
    SpecData->I_Mag_2D_unpolarized      = (float*) malloc(len);
    SpecData->I_Mag_2D_polarized        = (float*) malloc(len);
    SpecData->I_NucMag_2D               = (float*) malloc(len);
    SpecData->I_Mag_2D_spin_flip        = (float*) malloc(len);
    SpecData->I_Mag_2D_chiral           = (float*) malloc(len);
    SpecData->I_Mag_2D_spin_flip_pm     = (float*) malloc(len);
    SpecData->I_Mag_2D_spin_flip_mp     = (float*) malloc(len);
    SpecData->I_Mag_2D_non_spin_flip_pp = (float*) malloc(len);
    SpecData->I_Mag_2D_non_spin_flip_mm = (float*) malloc(len);
    SpecData->I_Mag_2D_sanspol_p        = (float*) malloc(len);
    SpecData->I_Mag_2D_sanspol_m        = (float*) malloc(len);

    // --- Initialize to zero ---
    memset(SpecData->I_Nuc_2D_unpolarized, 0, len);
    memset(SpecData->I_Mag_2D_unpolarized, 0, len);
    memset(SpecData->I_Mag_2D_polarized, 0, len);
    memset(SpecData->I_NucMag_2D, 0, len);
    memset(SpecData->I_Mag_2D_spin_flip, 0, len);
    memset(SpecData->I_Mag_2D_chiral, 0, len);
    memset(SpecData->I_Mag_2D_spin_flip_pm, 0, len);
    memset(SpecData->I_Mag_2D_spin_flip_mp, 0, len);
    memset(SpecData->I_Mag_2D_non_spin_flip_pp, 0, len);
    memset(SpecData->I_Mag_2D_non_spin_flip_mm, 0, len);
    memset(SpecData->I_Mag_2D_sanspol_p, 0, len);
    memset(SpecData->I_Mag_2D_sanspol_m, 0, len);
 
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

void free_SpectralData(AngularSpectrum* A)
{
    
}



