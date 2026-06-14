
.. _program_listing_file_src_gpu_kernels_AtomisticNucSANS_DiluteKernel.h:

Program Listing for File AtomisticNucSANS_DiluteKernel.h
========================================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_gpu_kernels_AtomisticNucSANS_DiluteKernel.h>` (``src/gpu_kernels/AtomisticNucSANS_DiluteKernel.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   
   // Computes purely nuclear atomistic SANS for a dilute system.
   // Active data layers:
   //   - NuclearData: local nuclear scattering length density data.
   // Inactive data layers:
   //   - MagnetizationData
   //   - StructureData
   //   - RotationData
   // Use case:
   //   Single-object or dilute-ensemble nuclear scattering without explicit
   //   object centers and without object-wise rotations.
   //
   // Placeholder for:
   // Atomistic_NucSANS_Kernel_dilute
   
   __global__ void Atomistic_NucSANS_Kernel_dilute(NuclearData NucData, ScatteringData SANSData) {
   
       // Input information:
       // N     : number of atoms
       // L     : number of points in Fourier space L = N_q*N_theta
       // K     : number of particles
       // x     : x-real-space position data in units of nano-meters
       // y     : y-real-space position data in units of nano-meters
       // z     : z-real-space position data in units of nano-meters
       // qy    : qy-Fourier-space coordinate in units of inverse nano-meters
       // qz    : qz-Fourier-space coordinate in units of inverse nano-meters
       // theta : theta-angle in Fourier space [theta = arctan2(qz, qy)] in radiant
   
       // output information:
   
       unsigned long int L = (*SANSData.N_q) * (*SANSData.N_theta);
       unsigned long int N_avg = *NucData.N_avg;
       unsigned long int W = *NucData.TotalAtomNumber;
   
       // float v = 1.0/((float)  (*NucData.K)) * powf(1.0/((float) (*NucData.N)), 2); // pre factor
       float v = 1.0 / ((float)W) * 1.0 / ((float)N_avg);
       ;
   
       int i = blockIdx.x * blockDim.x + threadIdx.x;
   
       float Nuc_real = 0.0;
       float Nuc_imag = 0.0;
       float Y = 0.0;
       float Z = 0.0;
   
       float Psi = 0.0;
   
       float cos_val = 0.0;
       float sin_val = 0.0;
   
       unsigned long int N_cum = 0;
       unsigned long int N_act = 0;
   
       if (i < L) {
           for (int k = 0; k < (*NucData.K); k++) {
   
               Nuc_real = 0.0;
               Nuc_imag = 0.0;
   
               N_cum = NucData.N_cum[k];
               N_act = NucData.N_act[k];
   
               for (int l = 0; l < N_act; l++) {
                   // atomic position composition
                   // X = MagData.RotMat[0] * MagData.x[l+k*N] \
                   //  + MagData.RotMat[3] * MagData.y[l+k*N] \
                   // + MagData.RotMat[6] * MagData.z[l+k*N];
                   Y = NucData.RotMat[1] * NucData.x[l + N_cum] + NucData.RotMat[4] * NucData.y[l + N_cum] +
                       NucData.RotMat[7] * NucData.z[l + N_cum];
                   Z = NucData.RotMat[2] * NucData.x[l + N_cum] + NucData.RotMat[5] * NucData.y[l + N_cum] +
                       NucData.RotMat[8] * NucData.z[l + N_cum];
   
                   // phase function
                   Psi = Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i];
   
                   // cosine and sine values
                   cos_val = cosf(Psi);
                   sin_val = sinf(Psi);
   
                   // cosine- and sine-summations
                   Nuc_real += NucData.Nuc[l + N_cum] * cos_val;
                   Nuc_imag -= NucData.Nuc[l + N_cum] * sin_val;
               }
   
               // nuclear SANS cross section projected in (qz, qy)-plane
               SANSData.S_Nuc_2D_unpolarized[i] += v * (Nuc_real * Nuc_real + Nuc_imag * Nuc_imag);
           }
       }
   }
