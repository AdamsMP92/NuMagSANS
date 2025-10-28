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
	// the total array length is 2 * Nq * (k_max + 1)
	float *I_Nuc_unpolarized;		  // nuclear SANS cross section
	float *I_Mag_unpolarized;		  // unpolarized magnetic SANS cross section
	float *I_Mag_polarized;			  // polarized magnetic SANS cross section
	float *I_NucMag;				  // nuclear-magnetic interference SANS cross section
	float *I_Mag_spin_flip;			  // spin-flip magnetic SANS cross section
	float *I_Mag_chiral;			  // chiral magnetic SANS cross section
	float *I_Mag_spin_flip_pm;		  // pm-spin-flip magnetic SANS cross section
	float *I_Mag_spin_flip_mp;		  // mp-spin-flip magnetic SANS cross section
	float *I_Mag_non_spin_flip_pp;	  // pp-non-spin-flip magnetic SANS cross section
	float *I_Mag_non_spin_flip_mm;	  // mm-non-spin-flip magnetic SANS cross section
	float *I_Mag_sanspol_p;			  // p-sanspol magnetic SANS cross section
	float *I_Mag_sanspol_m;			  // m-sanspol magnetic SANS cross section
};


void scale_SpectralData(ScalingFactors* ScalFactors, \
						SpectralData* SpecData, \
                        InputFileData* InputData){

	unsigned int Nq = *SpecData->Nq;
	unsigned int k_max = *SpecData->k_max;
	unsigned int L = 2 * Nq * (k_max + 1);

	for(unsigned long int l=0; l < L; l++){

		// basic SANS cross sections
		SpecData->I_Nuc_unpolarized[l] = SpecData->I_Nuc_unpolarized[l] * ScalFactors->Nuc_SANS_SF;
		SpecData->I_Mag_unpolarized[l] = SpecData->I_Mag_unpolarized[l] * ScalFactors->Mag_SANS_SF;
		SpecData->I_Mag_polarized[l] = SpecData->I_Mag_polarized[l] * ScalFactors->Mag_SANS_SF;;
		SpecData->I_NucMag[l] = SpecData->I_NucMag[l] * ScalFactors->NucMag_SANS_SF;
		SpecData->I_Mag_chiral[l] = SpecData->I_Mag_chiral[l] * ScalFactors->Mag_SANS_SF;

		// composed SANS cross sections
		SpecData->I_Mag_spin_flip[l] = SpecData->I_Mag_unpolarized[l] - SpecData->I_Mag_polarized[l];
		SpecData->I_Mag_spin_flip_pm[l] = SpecData->I_Mag_spin_flip[l] + SpecData->I_Mag_chiral[l];
		SpecData->I_Mag_spin_flip_mp[l] = SpecData->I_Mag_spin_flip[l] - SpecData->I_Mag_chiral[l];
		SpecData->I_Mag_non_spin_flip_pp[l] = SpecData->I_Nuc_unpolarized[l] + SpecData->I_NucMag[l] + SpecData->I_Mag_polarized[l];
		SpecData->I_Mag_non_spin_flip_mm[l] = SpecData->I_Nuc_unpolarized[l] - SpecData->I_NucMag[l] + SpecData->I_Mag_polarized[l];
		SpecData->I_Mag_sanspol_p[l] = SpecData->I_Mag_non_spin_flip_pp[l] + SpecData->I_Mag_spin_flip_pm[l];
		SpecData->I_Mag_sanspol_m[l] = SpecData->I_Mag_non_spin_flip_mm[l] + SpecData->I_Mag_spin_flip_mp[l];

	}
}

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

	cudaDeviceSynchronize();
	
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


void write2CSV_SpectralData(InputFileData* InputData,
                            SpectralData* SpecData,
                            int MagData_File_Index){
    LogSystem::write("");
    LogSystem::write("write spectral decomposition data to csv-files...");

    unsigned int Nq    = *SpecData->Nq;
    unsigned int k_max = *SpecData->k_max;
    unsigned int L     = (k_max + 1);
    unsigned int Ltot  = 2 * Nq * L;  // sin + cos arrays combined

    // === Ordner anlegen ===
    std::string foldername = InputData->SANSDataFoldername +
                             "/SANS_" + std::to_string(MagData_File_Index) + "/AngularSpectrum/";	
    mkdir(foldername.c_str(), 0777);

    auto write_component = [&](const std::string& fname, float* data) {
        std::ofstream fout(foldername + fname);
        if (!fout.is_open()) {
            LogSystem::write("Error opening " + fname);
            return;
        }

        // --- Header ---
        fout << "q";
        for (unsigned int k = 0; k <= k_max; ++k) fout << ",Ic_" << k;
        for (unsigned int k = 0; k <= k_max; ++k) fout << ",Is_" << k;
        fout << "\n";

        // --- Daten schreiben ---
        unsigned int offset_cos = 0;
        unsigned int offset_sin = Nq * L; // sin arrays beginnen nach allen cos Werten

        for (unsigned int i = 0; i < Nq; ++i) {
            fout << SpecData->q[i];

            // cos terms
            for (unsigned int k = 0; k <= k_max; ++k) {
                unsigned int idx = i + k * Nq + offset_cos;
                fout << "," << data[idx];
            }

            // sin terms
            for (unsigned int k = 0; k <= k_max; ++k) {
                unsigned int idx = i + k * Nq + offset_sin;
                fout << "," << data[idx];
            }

            fout << "\n";
        }

        fout.close();
        LogSystem::write(fname + " written.");
    };

    // === Write each spectrum ===
    write_component("I_Nuc_unpolarized.csv",      SpecData->I_Nuc_unpolarized);
    write_component("I_Mag_unpolarized.csv",      SpecData->I_Mag_unpolarized);
    write_component("I_Mag_polarized.csv",        SpecData->I_Mag_polarized);
    write_component("I_NucMag.csv",               SpecData->I_NucMag);
    write_component("I_Mag_chiral.csv",           SpecData->I_Mag_chiral);
    write_component("I_Mag_spin_flip.csv",        SpecData->I_Mag_spin_flip);
    write_component("I_Mag_spin_flip_pm.csv",     SpecData->I_Mag_spin_flip_pm);
    write_component("I_Mag_spin_flip_mp.csv",     SpecData->I_Mag_spin_flip_mp);
    write_component("I_Mag_non_spin_flip_pp.csv", SpecData->I_Mag_non_spin_flip_pp);
    write_component("I_Mag_non_spin_flip_mm.csv", SpecData->I_Mag_non_spin_flip_mm);
    write_component("I_Mag_sanspol_p.csv",        SpecData->I_Mag_sanspol_p);
    write_component("I_Mag_sanspol_m.csv",        SpecData->I_Mag_sanspol_m);

    LogSystem::write("All spectral CSV files written successfully.");
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



