
.. _program_listing_file_src_NuMagSANSlib_SpectralData.h:

Program Listing for File NuMagSANSlib_SpectralData.h
====================================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_NuMagSANSlib_SpectralData.h>` (``src/NuMagSANSlib_SpectralData.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

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
       float *I_Nuc_unpolarized;         // nuclear SANS cross section
       float *I_Mag_unpolarized;         // unpolarized magnetic SANS cross section
       float *I_Mag_polarized;           // polarized magnetic SANS cross section
       float *I_NucMag;                  // nuclear-magnetic interference SANS cross section
       float *I_Mag_spin_flip;           // spin-flip magnetic SANS cross section
       float *I_Mag_chiral;              // chiral magnetic SANS cross section
       float *I_Mag_spin_flip_pm;        // pm-spin-flip magnetic SANS cross section
       float *I_Mag_spin_flip_mp;        // mp-spin-flip magnetic SANS cross section
       float *I_Mag_non_spin_flip_pp;    // pp-non-spin-flip magnetic SANS cross section
       float *I_Mag_non_spin_flip_mm;    // mm-non-spin-flip magnetic SANS cross section
       float *I_Mag_sanspol_p;           // p-sanspol magnetic SANS cross section
       float *I_Mag_sanspol_m;           // m-sanspol magnetic SANS cross section
   
       // normalized angular amplitude spectrum
       // the total array length is 2 * (k_max + 1)
       float *A_Nuc_unpolarized;         // nuclear SANS cross section
       float *A_Mag_unpolarized;         // unpolarized magnetic SANS cross section
       float *A_Mag_polarized;           // polarized magnetic SANS cross section
       float *A_NucMag;                  // nuclear-magnetic interference SANS cross section
       float *A_Mag_spin_flip;           // spin-flip magnetic SANS cross section
       float *A_Mag_chiral;              // chiral magnetic SANS cross section
       float *A_Mag_spin_flip_pm;        // pm-spin-flip magnetic SANS cross section
       float *A_Mag_spin_flip_mp;        // mp-spin-flip magnetic SANS cross section
       float *A_Mag_non_spin_flip_pp;    // pp-non-spin-flip magnetic SANS cross section
       float *A_Mag_non_spin_flip_mm;    // mm-non-spin-flip magnetic SANS cross section
       float *A_Mag_sanspol_p;           // p-sanspol magnetic SANS cross section
       float *A_Mag_sanspol_m;           // m-sanspol magnetic SANS cross section
   
   };
   
   
   void scale_SpectralData(ScalingFactors* ScalFactors, SpectralData* SpecData, InputFileData* InputData){
   
       unsigned int Nq = *SpecData->Nq;
       unsigned int k_max = *SpecData->k_max;
       unsigned int L = 2 * Nq * (k_max + 1);
       unsigned int Lk = 2 * (k_max + 1);
   
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
   
       // =====================================================
       // (2) scaling and normalization of the amplitude spectrum
       // =====================================================
       auto normalizeAmplitude = [&](float* A_array, float scaleFactor) {
           // (a) Scaling
           for (unsigned int k = 0; k < Lk; ++k)
               A_array[k] *= scaleFactor;
   
           // (b) compute norm of the ampliutde spectrum
           float normSum = 0.0f;
           for (unsigned int k = 0; k <= k_max; ++k) {
               float Ac = A_array[k];
               float As = A_array[k + (k_max + 1)];
               normSum += sqrtf(Ac * Ac + As * As);
           }
   
           if (normSum == 0.0f) return;
   
           // (c) normalization
           for (unsigned int k = 0; k <= k_max; ++k) {
               A_array[k]               /= normSum;
               A_array[k + (k_max + 1)] /= normSum;
           }
       };
   
       // normalizing of the base amplitudes
       normalizeAmplitude(SpecData->A_Nuc_unpolarized, ScalFactors->Nuc_SANS_SF);
       normalizeAmplitude(SpecData->A_Mag_unpolarized, ScalFactors->Mag_SANS_SF);
       normalizeAmplitude(SpecData->A_Mag_polarized,   ScalFactors->Mag_SANS_SF);
       normalizeAmplitude(SpecData->A_NucMag,          ScalFactors->NucMag_SANS_SF);
       normalizeAmplitude(SpecData->A_Mag_chiral,      ScalFactors->Mag_SANS_SF);
   
       // =====================================================
       // composition of the secondary amplitudes
       // =====================================================
       for (unsigned int k = 0; k < Lk; ++k) {
           SpecData->A_Mag_spin_flip[k]        = SpecData->A_Mag_unpolarized[k] - SpecData->A_Mag_polarized[k];
           SpecData->A_Mag_spin_flip_pm[k]     = SpecData->A_Mag_spin_flip[k] + SpecData->A_Mag_chiral[k];
           SpecData->A_Mag_spin_flip_mp[k]     = SpecData->A_Mag_spin_flip[k] - SpecData->A_Mag_chiral[k];
           SpecData->A_Mag_non_spin_flip_pp[k] = SpecData->A_Nuc_unpolarized[k] + SpecData->A_NucMag[k] + SpecData->A_Mag_polarized[k];
           SpecData->A_Mag_non_spin_flip_mm[k] = SpecData->A_Nuc_unpolarized[k] - SpecData->A_NucMag[k] + SpecData->A_Mag_polarized[k];
           SpecData->A_Mag_sanspol_p[k]        = SpecData->A_Mag_non_spin_flip_pp[k] + SpecData->A_Mag_spin_flip_pm[k];
           SpecData->A_Mag_sanspol_m[k]        = SpecData->A_Mag_non_spin_flip_mm[k] + SpecData->A_Mag_spin_flip_mp[k];
       }
   
       // normalize the secondary amplitudes
       normalizeAmplitude(SpecData->A_Mag_spin_flip,        1.0f);
       normalizeAmplitude(SpecData->A_Mag_spin_flip_pm,     1.0f);
       normalizeAmplitude(SpecData->A_Mag_spin_flip_mp,     1.0f);
       normalizeAmplitude(SpecData->A_Mag_non_spin_flip_pp, 1.0f);
       normalizeAmplitude(SpecData->A_Mag_non_spin_flip_mm, 1.0f);
       normalizeAmplitude(SpecData->A_Mag_sanspol_p,        1.0f);
       normalizeAmplitude(SpecData->A_Mag_sanspol_m,        1.0f);
   
   }
   
   void allocate_SpectralData_RAM(InputFileData* InputData, ScatteringData* SANSData, SpectralData* SpecData){
   
       SpecData->Nq = (unsigned int*) malloc(sizeof(unsigned int));
       SpecData->Ntheta = (unsigned int*) malloc(sizeof(unsigned int));
       SpecData->k_max = (unsigned int*) malloc(sizeof(unsigned int));
       SpecData->dtheta = (float*) malloc(sizeof(float));
    
       *SpecData->Nq = InputData->N_q;
       *SpecData->Ntheta = InputData->N_theta;
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
   
        // allocate spectral arrays
       size_t lenA = 2 *  ((*(SpecData->k_max)) + 1) * sizeof(float);
       SpecData->A_Nuc_unpolarized      = (float*) malloc(lenA);
       SpecData->A_Mag_unpolarized      = (float*) malloc(lenA);
       SpecData->A_Mag_polarized        = (float*) malloc(lenA);
       SpecData->A_NucMag               = (float*) malloc(lenA);
       SpecData->A_Mag_spin_flip        = (float*) malloc(lenA);
       SpecData->A_Mag_chiral           = (float*) malloc(lenA);
       SpecData->A_Mag_spin_flip_pm     = (float*) malloc(lenA);
       SpecData->A_Mag_spin_flip_mp     = (float*) malloc(lenA);
       SpecData->A_Mag_non_spin_flip_pp = (float*) malloc(lenA);
       SpecData->A_Mag_non_spin_flip_mm = (float*) malloc(lenA);
       SpecData->A_Mag_sanspol_p        = (float*) malloc(lenA);
       SpecData->A_Mag_sanspol_m        = (float*) malloc(lenA);
   
       // --- Initialize to zero ---
       memset(SpecData->A_Nuc_unpolarized, 0, lenA);
       memset(SpecData->A_Mag_unpolarized, 0, lenA);
       memset(SpecData->A_Mag_polarized, 0, lenA);
       memset(SpecData->A_NucMag, 0, lenA);
       memset(SpecData->A_Mag_spin_flip, 0, lenA);
       memset(SpecData->A_Mag_chiral, 0, lenA);
       memset(SpecData->A_Mag_spin_flip_pm, 0, lenA);
       memset(SpecData->A_Mag_spin_flip_mp, 0, lenA);
       memset(SpecData->A_Mag_non_spin_flip_pp, 0, lenA);
       memset(SpecData->A_Mag_non_spin_flip_mm, 0, lenA);
       memset(SpecData->A_Mag_sanspol_p, 0, lenA);
       memset(SpecData->A_Mag_sanspol_m, 0, lenA);
    
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
       size_t len_A = 2 * (k_max + 1) * sizeof(float);
   
       // Allocate arrays on GPU
       cudaMalloc(&(SpecData_gpu->q), len_q);
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
   
   
       cudaMalloc(&(SpecData_gpu->A_Nuc_unpolarized), len_A);
       cudaMalloc(&(SpecData_gpu->A_Mag_unpolarized), len_A);
       cudaMalloc(&(SpecData_gpu->A_Mag_polarized), len_A);
       cudaMalloc(&(SpecData_gpu->A_NucMag), len_A);
       cudaMalloc(&(SpecData_gpu->A_Mag_spin_flip), len_A);
       cudaMalloc(&(SpecData_gpu->A_Mag_chiral), len_A);
       cudaMalloc(&(SpecData_gpu->A_Mag_spin_flip_pm), len_A);
       cudaMalloc(&(SpecData_gpu->A_Mag_spin_flip_mp), len_A);
       cudaMalloc(&(SpecData_gpu->A_Mag_non_spin_flip_pp), len_A);
       cudaMalloc(&(SpecData_gpu->A_Mag_non_spin_flip_mm), len_A);
       cudaMalloc(&(SpecData_gpu->A_Mag_sanspol_p), len_A);
       cudaMalloc(&(SpecData_gpu->A_Mag_sanspol_m), len_A);
   
       // Initialize GPU arrays with zeros
       cudaMemset(SpecData_gpu->A_Nuc_unpolarized,      0, len_A);
       cudaMemset(SpecData_gpu->A_Mag_unpolarized,      0, len_A);
       cudaMemset(SpecData_gpu->A_Mag_polarized,        0, len_A);
       cudaMemset(SpecData_gpu->A_NucMag,               0, len_A);
       cudaMemset(SpecData_gpu->A_Mag_spin_flip,        0, len_A);
       cudaMemset(SpecData_gpu->A_Mag_chiral,           0, len_A);
       cudaMemset(SpecData_gpu->A_Mag_spin_flip_pm,     0, len_A);
       cudaMemset(SpecData_gpu->A_Mag_spin_flip_mp,     0, len_A);
       cudaMemset(SpecData_gpu->A_Mag_non_spin_flip_pp, 0, len_A);
       cudaMemset(SpecData_gpu->A_Mag_non_spin_flip_mm, 0, len_A);
       cudaMemset(SpecData_gpu->A_Mag_sanspol_p,        0, len_A);
       cudaMemset(SpecData_gpu->A_Mag_sanspol_m,        0, len_A);
   
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
       size_t len_A = 2 * (k_max + 1) * sizeof(float);
   
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
   
       // Copy spectral arrays 
       cudaMemcpy(S->A_Nuc_unpolarized,      S_gpu->A_Nuc_unpolarized,      len_A, cudaMemcpyDeviceToHost);
       cudaMemcpy(S->A_Mag_unpolarized,      S_gpu->A_Mag_unpolarized,      len_A, cudaMemcpyDeviceToHost);
       cudaMemcpy(S->A_Mag_polarized,        S_gpu->A_Mag_polarized,        len_A, cudaMemcpyDeviceToHost);
       cudaMemcpy(S->A_NucMag,               S_gpu->A_NucMag,               len_A, cudaMemcpyDeviceToHost);
       cudaMemcpy(S->A_Mag_spin_flip,        S_gpu->A_Mag_spin_flip,        len_A, cudaMemcpyDeviceToHost);
       cudaMemcpy(S->A_Mag_chiral,           S_gpu->A_Mag_chiral,           len_A, cudaMemcpyDeviceToHost);
       cudaMemcpy(S->A_Mag_spin_flip_pm,     S_gpu->A_Mag_spin_flip_pm,     len_A, cudaMemcpyDeviceToHost);
       cudaMemcpy(S->A_Mag_spin_flip_mp,     S_gpu->A_Mag_spin_flip_mp,     len_A, cudaMemcpyDeviceToHost);
       cudaMemcpy(S->A_Mag_non_spin_flip_pp, S_gpu->A_Mag_non_spin_flip_pp, len_A, cudaMemcpyDeviceToHost);
       cudaMemcpy(S->A_Mag_non_spin_flip_mm, S_gpu->A_Mag_non_spin_flip_mm, len_A, cudaMemcpyDeviceToHost);
       cudaMemcpy(S->A_Mag_sanspol_p,        S_gpu->A_Mag_sanspol_p,        len_A, cudaMemcpyDeviceToHost);
       cudaMemcpy(S->A_Mag_sanspol_m,        S_gpu->A_Mag_sanspol_m,        len_A, cudaMemcpyDeviceToHost);
   
       // === 4. Synchronize to ensure data is ready ===
       cudaDeviceSynchronize();
   
   }
   
   
   // void write2CSV_SpectralData(InputFileData* InputData, \
   //                             SpectralData* SpecData, \
   //                             int MagData_File_Index){
   //     LogSystem::write("");
   //     LogSystem::write("write spectral decomposition data to csv-files...");
   
   //     unsigned int Nq    = *SpecData->Nq;
   //     unsigned int k_max = *SpecData->k_max;
   //     unsigned int L     = (k_max + 1);
   //     //unsigned int Ltot  = 2 * Nq * L;  // sin + cos arrays combined
   
   //     // === Ordner anlegen ===
   //     std::string foldername = InputData->SANSDataFoldername +
   //                              "/SANS_" + std::to_string(MagData_File_Index) + "/AngularSpectrum/";    
   //     mkdir(foldername.c_str(), 0777);
   
   //  // =====================================================
   //     // (1) Hilfsfunktion: Spektral-Intensitäten I(q,k)
   //     // =====================================================
   
   //     auto write_component = [&](const std::string& fname, float* data) {
   //         std::ofstream fout(foldername + fname);
   //         if (!fout.is_open()) {
   //             LogSystem::write("Error opening " + fname);
   //             return;
   //         }
   
   //         // --- Header ---
   //         fout << "q";
   //         for (unsigned int k = 0; k <= k_max; ++k) fout << ",Ic_" << k;
   //         for (unsigned int k = 0; k <= k_max; ++k) fout << ",Is_" << k;
   //         fout << "\n";
   
   //         // --- Daten schreiben ---
   //         unsigned int offset_cos = 0;
   //         unsigned int offset_sin = Nq * L; // sin arrays beginnen nach allen cos Werten
   
   //         for (unsigned int i = 0; i < Nq; ++i) {
   //             fout << SpecData->q[i];
   
   //             // cos terms
   //             for (unsigned int k = 0; k <= k_max; ++k) {
   //                 unsigned int idx = i + k * Nq + offset_cos;
   //                 fout << "," << data[idx];
   //             }
   
   //             // sin terms
   //             for (unsigned int k = 0; k <= k_max; ++k) {
   //                 unsigned int idx = i + k * Nq + offset_sin;
   //                 fout << "," << data[idx];
   //             }
   
   //             fout << "\n";
   //         }
   
   //         fout.close();
   //         LogSystem::write(fname + " written.");
   //     };
   
   
   //  // =====================================================
   //     // (2) Hilfsfunktion: Amplituden A(k)
   //     // =====================================================
   //     auto write_amplitude = [&](const std::string& fname, float* A_array) {
   //         std::ofstream fout(foldername + fname);
   //         if (!fout.is_open()) {
   //             LogSystem::write("Error opening " + fname);
   //             return;
   //         }
   
   //         // --- Header ---
   //         fout << "k,A_cos,A_sin\n";
   
   //         for (unsigned int k = 0; k <= k_max; ++k) {
   //             float A_cos = A_array[k];
   //             float A_sin = A_array[k + (k_max + 1)];
   //             fout << k << "," << A_cos << "," << A_sin << "\n";
   //         }
   
   //         fout.close();
   //         LogSystem::write(fname + " written.");
   //     };
   
   //     std::string intensityFolder = foldername + "Intensities/";
   //     mkdir(intensityFolder.c_str(), 0777);
   
   //     // === Write each spectrum ===
   //     write_component("Intensities/I_Nuc_unpolarized.csv",      SpecData->I_Nuc_unpolarized);
   //     write_component("Intensities/I_Mag_unpolarized.csv",      SpecData->I_Mag_unpolarized);
   //     write_component("Intensities/I_Mag_polarized.csv",        SpecData->I_Mag_polarized);
   //     write_component("Intensities/I_NucMag.csv",               SpecData->I_NucMag);
   //     write_component("Intensities/I_Mag_chiral.csv",           SpecData->I_Mag_chiral);
   //     write_component("Intensities/I_Mag_spin_flip.csv",        SpecData->I_Mag_spin_flip);
   //     write_component("Intensities/I_Mag_spin_flip_pm.csv",     SpecData->I_Mag_spin_flip_pm);
   //     write_component("Intensities/I_Mag_spin_flip_mp.csv",     SpecData->I_Mag_spin_flip_mp);
   //     write_component("Intensities/I_Mag_non_spin_flip_pp.csv", SpecData->I_Mag_non_spin_flip_pp);
   //     write_component("Intensities/I_Mag_non_spin_flip_mm.csv", SpecData->I_Mag_non_spin_flip_mm);
   //     write_component("Intensities/I_Mag_sanspol_p.csv",        SpecData->I_Mag_sanspol_p);
   //     write_component("Intensities/I_Mag_sanspol_m.csv",        SpecData->I_Mag_sanspol_m);
   
   
   //  // =====================================================
   //     // (4) Export der Angular-Amplituden A(k)
   //     // =====================================================
   //     std::string ampFolder = foldername + "Amplitudes/";
   //     mkdir(ampFolder.c_str(), 0777);
   
   //     write_amplitude("Amplitudes/A_Nuc_unpolarized.csv",      SpecData->A_Nuc_unpolarized);
   //     write_amplitude("Amplitudes/A_Mag_unpolarized.csv",      SpecData->A_Mag_unpolarized);
   //     write_amplitude("Amplitudes/A_Mag_polarized.csv",        SpecData->A_Mag_polarized);
   //     write_amplitude("Amplitudes/A_NucMag.csv",               SpecData->A_NucMag);
   //     write_amplitude("Amplitudes/A_Mag_chiral.csv",           SpecData->A_Mag_chiral);
   //     write_amplitude("Amplitudes/A_Mag_spin_flip.csv",        SpecData->A_Mag_spin_flip);
   //     write_amplitude("Amplitudes/A_Mag_spin_flip_pm.csv",     SpecData->A_Mag_spin_flip_pm);
   //     write_amplitude("Amplitudes/A_Mag_spin_flip_mp.csv",     SpecData->A_Mag_spin_flip_mp);
   //     write_amplitude("Amplitudes/A_Mag_non_spin_flip_pp.csv", SpecData->A_Mag_non_spin_flip_pp);
   //     write_amplitude("Amplitudes/A_Mag_non_spin_flip_mm.csv", SpecData->A_Mag_non_spin_flip_mm);
   //     write_amplitude("Amplitudes/A_Mag_sanspol_p.csv",        SpecData->A_Mag_sanspol_p);
   //     write_amplitude("Amplitudes/A_Mag_sanspol_m.csv",        SpecData->A_Mag_sanspol_m);
   
   //     LogSystem::write("All spectral CSV files written successfully.");
   // }
   
   
   // #################################################################################################
   // #################################################################################################
   
   struct SpectralComponent {
       const char* name;
       const float* data;   // raw pointer
   };
   
   
   void write_spectral_csv(
       const std::string& folder,
       const SpectralComponent& comp,
       SpectralData* SpecData)
   {
       unsigned int Nq    = *SpecData->Nq;
       unsigned int k_max = *SpecData->k_max;
       unsigned int L     = k_max + 1;
   
       std::ofstream fout(folder + std::string(comp.name) + ".csv");
   
       // Header
       fout << "q";
       for (unsigned int k = 0; k <= k_max; ++k) fout << ",Ic_" << k;
       for (unsigned int k = 0; k <= k_max; ++k) fout << ",Is_" << k;
       fout << "\n";
   
       unsigned int offset_cos = 0;
       unsigned int offset_sin = Nq * L;
   
       for (unsigned int i = 0; i < Nq; ++i)
       {
           fout << SpecData->q[i];
   
           for (unsigned int k = 0; k <= k_max; ++k)
               fout << "," << comp.data[i + k * Nq + offset_cos];
   
           for (unsigned int k = 0; k <= k_max; ++k)
               fout << "," << comp.data[i + k * Nq + offset_sin];
   
           fout << "\n";
       }
   }
   
   
   void write_amplitude_csv(
       const std::string& folder,
       const SpectralComponent& comp,
       SpectralData* SpecData)
   {
       unsigned int k_max = *SpecData->k_max;
       unsigned int L     = k_max + 1;
   
       std::ofstream fout(folder + std::string(comp.name) + ".csv");
       if(!fout.is_open()){
           LogSystem::write("Error opening amplitude file: " + std::string(comp.name));
           return;
       }
   
       fout << "k,A_cos,A_sin\n";
   
       for(unsigned int k = 0; k <= k_max; ++k)
       {
           float A_cos = comp.data[k];
           float A_sin = comp.data[k + L];
           fout << k << "," << A_cos << "," << A_sin << "\n";
       }
   }
   
   
   std::vector<SpectralComponent> build_spectral_intensities(SpectralData* SpecData)
   {
       return {
           {"I_Nuc_unpolarized",      SpecData->I_Nuc_unpolarized},
           {"I_Mag_unpolarized",      SpecData->I_Mag_unpolarized},
           {"I_Mag_polarized",        SpecData->I_Mag_polarized},
           {"I_NucMag",               SpecData->I_NucMag},
           {"I_Mag_chiral",           SpecData->I_Mag_chiral},
           {"I_Mag_spin_flip",        SpecData->I_Mag_spin_flip},
           {"I_Mag_spin_flip_pm",     SpecData->I_Mag_spin_flip_pm},
           {"I_Mag_spin_flip_mp",     SpecData->I_Mag_spin_flip_mp},
           {"I_Mag_non_spin_flip_pp", SpecData->I_Mag_non_spin_flip_pp},
           {"I_Mag_non_spin_flip_mm", SpecData->I_Mag_non_spin_flip_mm},
           {"I_Mag_sanspol_p",        SpecData->I_Mag_sanspol_p},
           {"I_Mag_sanspol_m",        SpecData->I_Mag_sanspol_m}
       };
   }
   
   std::vector<SpectralComponent> build_spectral_amplitudes(SpectralData* SpecData)
   {
       return {
           {"A_Nuc_unpolarized",      SpecData->A_Nuc_unpolarized},
           {"A_Mag_unpolarized",      SpecData->A_Mag_unpolarized},
           {"A_Mag_polarized",        SpecData->A_Mag_polarized},
           {"A_NucMag",               SpecData->A_NucMag},
           {"A_Mag_chiral",           SpecData->A_Mag_chiral},
           {"A_Mag_spin_flip",        SpecData->A_Mag_spin_flip},
           {"A_Mag_spin_flip_pm",     SpecData->A_Mag_spin_flip_pm},
           {"A_Mag_spin_flip_mp",     SpecData->A_Mag_spin_flip_mp},
           {"A_Mag_non_spin_flip_pp", SpecData->A_Mag_non_spin_flip_pp},
           {"A_Mag_non_spin_flip_mm", SpecData->A_Mag_non_spin_flip_mm},
           {"A_Mag_sanspol_p",        SpecData->A_Mag_sanspol_p},
           {"A_Mag_sanspol_m",        SpecData->A_Mag_sanspol_m}
       };
   }
   
   
   
   void write2CSV_SpectralData(
       InputFileData* InputData,
       SpectralData* SpecData,
       int MagData_File_Index)
   {
       LogSystem::write("");
       LogSystem::write("write spectral decomposition data to csv-files...");
   
       std::string baseFolder =
           InputData->SANSDataFoldername +
           "/SANS_" + std::to_string(MagData_File_Index) +
           "/AngularSpectrum/";
   
       mkdir(baseFolder.c_str(), 0777);
   
       // ================================
       // Intensities
       // ================================
       std::string intensityFolder = baseFolder + "Intensities/";
       mkdir(intensityFolder.c_str(), 0777);
   
       auto intensities = build_spectral_intensities(SpecData);
   
       for(const auto& comp : intensities)
           write_spectral_csv(intensityFolder, comp, SpecData);
   
       // ================================
       // Amplitudes
       // ================================
       std::string ampFolder = baseFolder + "Amplitudes/";
       mkdir(ampFolder.c_str(), 0777);
   
       auto amplitudes = build_spectral_amplitudes(SpecData);
   
       for(const auto& comp : amplitudes)
           write_amplitude_csv(ampFolder, comp, SpecData);
   
       LogSystem::write("All spectral CSV files written successfully.");
   }
   
   
   
   
   
   
   
   
   
   
   // #################################################################################################
   // #################################################################################################
   
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
       
       free(S->A_Nuc_unpolarized);
       free(S->A_Mag_unpolarized);
       free(S->A_Mag_polarized);
       free(S->A_NucMag);
       free(S->A_Mag_spin_flip);
       free(S->A_Mag_chiral);
       free(S->A_Mag_spin_flip_pm);
       free(S->A_Mag_spin_flip_mp);
       free(S->A_Mag_non_spin_flip_pp);
       free(S->A_Mag_non_spin_flip_mm);
       free(S->A_Mag_sanspol_p);
       free(S->A_Mag_sanspol_m);
   
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
   
       cudaFree(S_gpu->A_Nuc_unpolarized); 
       cudaDeviceSynchronize();
       cudaFree(S_gpu->A_Mag_unpolarized);
       cudaDeviceSynchronize();
       cudaFree(S_gpu->A_Mag_polarized); 
       cudaDeviceSynchronize();
       cudaFree(S_gpu->A_NucMag);
       cudaDeviceSynchronize();
       cudaFree(S_gpu->A_Mag_spin_flip); 
       cudaDeviceSynchronize();
       cudaFree(S_gpu->A_Mag_chiral);
       cudaDeviceSynchronize();
       cudaFree(S_gpu->A_Mag_spin_flip_pm); 
       cudaDeviceSynchronize();
       cudaFree(S_gpu->A_Mag_spin_flip_mp);
       cudaDeviceSynchronize();
       cudaFree(S_gpu->A_Mag_non_spin_flip_pp); 
       cudaDeviceSynchronize();
       cudaFree(S_gpu->A_Mag_non_spin_flip_mm);
       cudaDeviceSynchronize();
       cudaFree(S_gpu->A_Mag_sanspol_p); 
       cudaDeviceSynchronize();
       cudaFree(S_gpu->A_Mag_sanspol_m);
       cudaDeviceSynchronize();
   
   }
   
   
   
