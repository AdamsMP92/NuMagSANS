// File         : NuMagSANSlib_gpuKernel.h
// Author       : Michael Philipp ADAMS, M.Sc.
// Company      : University of Luxembourg
// Department   : Department of Physics and Materials Sciences
// Group        : NanoMagnetism Group
// Group Leader : Prof. Andreas Michels
// Version      : 28 November 2024
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





#include <cuda_runtime.h>
#include <math_constants.h>  // für M_PI

// ============================================================================
// GPU Kernel: compute angular spectra from 2D scattering data
// ============================================================================
__global__
void ComputeSpectralDecomposition(ScatteringData SANSData,
                                  SpectralData SpecData)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int Nq     = *SANSData.N_q;
    unsigned int Ntheta = *SANSData.N_theta;
    unsigned int k_max  = *SpecData.k_max;
    float dtheta        = *SANSData.dtheta;

    if (i >= Nq) return;

    // Loop over spectral modes k = 0...k_max-1
    for (unsigned int k = 0; k < k_max; ++k) {

        float sum_Nuc = 0.0f;
        float sum_Mag = 0.0f;
        float sum_Pol = 0.0f;
        float sum_NucMag = 0.0f;
        float sum_Chiral = 0.0f;

        // --- integrate over theta using trapezoidal rule ---
        for (unsigned int j = 0; j < Ntheta - 1; ++j) {
            float w = 0.5f * dtheta;  // Trapezregel-Gewicht
            float theta1 = j * dtheta;
            float theta2 = (j + 1) * dtheta;

            // Beispiel: spektrale Gewichtung (Fourier-Komponente)
            float phi_k1 = cosf(k * theta1);
            float phi_k2 = cosf(k * theta2);

            sum_Nuc    += w * (SANSData.S_Nuc_2D_unpolarized[j + i * Ntheta] * phi_k1 +
                               SANSData.S_Nuc_2D_unpolarized[j + 1 + i * Ntheta] * phi_k2);
            sum_Mag    += w * (SANSData.S_Mag_2D_unpolarized[j + i * Ntheta] * phi_k1 +
                               SANSData.S_Mag_2D_unpolarized[j + 1 + i * Ntheta] * phi_k2);
            sum_Pol    += w * (SANSData.S_Mag_2D_polarized[j + i * Ntheta] * phi_k1 +
                               SANSData.S_Mag_2D_polarized[j + 1 + i * Ntheta] * phi_k2);
            sum_NucMag += w * (SANSData.S_NucMag_2D[j + i * Ntheta] * phi_k1 +
                               SANSData.S_NucMag_2D[j + 1 + i * Ntheta] * phi_k2);
            sum_Chiral += w * (SANSData.S_Mag_2D_chiral[j + i * Ntheta] * phi_k1 +
                               SANSData.S_Mag_2D_chiral[j + 1 + i * Ntheta] * phi_k2);
        }

        // --- Ergebnis in Spektraldaten schreiben ---
        unsigned int idx = i * k_max + k;
        SpecData.I_Nuc_unpolarized[idx] = sum_Nuc / (4.0f * (float)M_PI);
        SpecData.I_Mag_unpolarized[idx] = sum_Mag / (4.0f * (float)M_PI);
        SpecData.I_Mag_polarized[idx]   = sum_Pol / (4.0f * (float)M_PI);
        SpecData.I_NucMag[idx]          = sum_NucMag / (4.0f * (float)M_PI);
        SpecData.I_Mag_chiral[idx]      = sum_Chiral / (4.0f * (float)M_PI);
    }
}




// GPU Kernel for the computation of the azimuthally averaged SANS cross section /////////////////////////////////////////////////////////////
// integration using trapezoidal rule ////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__
void AzimuthalAverage(ScatteringData SANSData){

	int i = blockIdx.x * blockDim.x + threadIdx.x;

	unsigned int N_theta = *SANSData.N_theta;
    unsigned int N_q = *SANSData.N_q;
    float dtheta = *SANSData.dtheta;
  
    if(i < N_q){
         for(int j=0; j<N_theta-1; j++){

             SANSData.S_Nuc_1D_unpolarized[i] += SANSData.S_Nuc_2D_unpolarized[j + i*N_theta] \
                                               + SANSData.S_Nuc_2D_unpolarized[j + i*N_theta + 1];
             SANSData.S_Mag_1D_unpolarized[i] += SANSData.S_Mag_2D_unpolarized[j + i*N_theta] \
                                               + SANSData.S_Mag_2D_unpolarized[j + i*N_theta + 1];
			 SANSData.S_Mag_1D_polarized[i] += SANSData.S_Mag_2D_polarized[j + i*N_theta] \
			  							     + SANSData.S_Mag_2D_polarized[j + i*N_theta + 1];
			 SANSData.S_NucMag_1D[i] += SANSData.S_NucMag_2D[j + i*N_theta] \
									  + SANSData.S_NucMag_2D[j + i*N_theta + 1];
			 SANSData.S_Mag_1D_chiral[i] += SANSData.S_Mag_2D_chiral[j + i*N_theta] \
			 						      + SANSData.S_Mag_2D_chiral[j + i*N_theta + 1];

          }
          
         SANSData.S_Nuc_1D_unpolarized[i] = SANSData.S_Nuc_1D_unpolarized[i]/(4.0*M_PI)*dtheta;
         SANSData.S_Mag_1D_unpolarized[i] = SANSData.S_Mag_1D_unpolarized[i]/(4.0*M_PI)*dtheta;
		 SANSData.S_Mag_1D_polarized[i] = SANSData.S_Mag_1D_polarized[i]/(4.0*M_PI)*dtheta;
		 SANSData.S_NucMag_1D[i] = SANSData.S_NucMag_1D[i]/(4.0*M_PI)*dtheta;
		 SANSData.S_Mag_1D_chiral[i] = SANSData.S_Mag_1D_chiral[i]/(4.0*M_PI)*dtheta;

     }
}



// computes in the first step the correlation function c(r) and the by multiplication with r^2 the pair-distance function
 // here we take into account the limit of sin(x)/x at x-> 0 and so the singularity is fixed
 __global__
 void DistributionFunctions(ScatteringData SANSData){


	// spherical hankel transform using a trapezoidal integration rule
  
       int i = blockIdx.x * blockDim.x + threadIdx.x;

       unsigned int N_r = *SANSData.N_r;
       unsigned int N_q = *SANSData.N_q;
       float dq = *SANSData.dq;
       float qr1 = 0.0;
       float qr2 = 0.0;
       bool b1 = false; 
       bool b2 = false; 
       float s1 = 0.0;
       float s2 = 0.0;
  
       if(i < N_r){
               for(int j=0; j<N_q-1; j++){
  
                   qr1 = SANSData.q_1D[j] * SANSData.r_1D[i];
                   b1 = (qr1 == 0.0f);
                   s1 = (sin(qr1)/(qr1 + (float)b1) + (float)b1) * pow(SANSData.q_1D[j], 2);
                   
	               qr2 = SANSData.q_1D[j+1] * SANSData.r_1D[i];
	               b2 = (qr2 == 0.0f);
        	       s2 = (sin(qr2)/(qr2 + (float)b2) + (float)b2) * pow(SANSData.q_1D[j+1], 2);

             	   SANSData.c_Nuc_unpolarized[i] += SANSData.S_Nuc_1D_unpolarized[j]  * s1 \
             	   								  + SANSData.S_Nuc_1D_unpolarized[j+1] * s2;
             	   								  
               	   SANSData.c_Mag_unpolarized[i] += SANSData.S_Mag_1D_unpolarized[j]  * s1 \
               	   								  + SANSData.S_Mag_1D_unpolarized[j+1] * s2;

            	   SANSData.c_NucMag[i] += SANSData.S_NucMag_1D[j]  * s1 \
         	               			     + SANSData.S_NucMag_1D[j+1] * s2;
               	   								  
               	   SANSData.c_Mag_polarized[i] += SANSData.S_Mag_1D_polarized[j]  * s1 \
               	   							    + SANSData.S_Mag_1D_polarized[j+1] * s2;

               	   SANSData.c_Mag_chiral[i] += SANSData.S_Mag_1D_chiral[j]  * s1 \
               	                  	   		 + SANSData.S_Mag_1D_chiral[j+1] * s2;
               	   							    
           		}

		        SANSData.c_Nuc_unpolarized[i] = SANSData.c_Nuc_unpolarized[i]/2.0 * dq;
                SANSData.p_Nuc_unpolarized[i] = SANSData.c_Nuc_unpolarized[i] * pow(SANSData.r_1D[i], 2);

                SANSData.c_Mag_unpolarized[i] = SANSData.c_Mag_unpolarized[i]/2.0 * dq;
                SANSData.p_Mag_unpolarized[i] = SANSData.c_Mag_unpolarized[i] * pow(SANSData.r_1D[i], 2);

                SANSData.c_NucMag[i] = SANSData.c_NucMag[i]/2.0 * dq;
                SANSData.p_NucMag[i] = SANSData.c_NucMag[i] * pow(SANSData.r_1D[i], 2);
   
                SANSData.c_Mag_polarized[i] = SANSData.c_Mag_polarized[i]/2.0 * dq;
                SANSData.p_Mag_polarized[i] = SANSData.c_Mag_polarized[i] * pow(SANSData.r_1D[i], 2);

                SANSData.c_Mag_chiral[i] = SANSData.c_Mag_chiral[i]/2.0 * dq;
                SANSData.p_Mag_chiral[i] = SANSData.c_Mag_chiral[i] * pow(SANSData.r_1D[i], 2);

      }
 }



  // 2D correlation functions #################################################################################
__global__
void CorrelationFunction_2D(ScatteringData SANSData){
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    int L_fourier = (*SANSData.N_q) * (*SANSData.N_theta);
    float v = 1.0/((float) L_fourier);
    int L_real = (*SANSData.N_r) * (*SANSData.N_alpha);
    float c = 0.0;

    if(i < L_real){
        for(int k = 0; k<L_fourier; k++){
            c = cos(SANSData.qy_2D[k] * SANSData.ry_2D[i] + SANSData.qz_2D[k] * SANSData.rz_2D[i]);
            SANSData.Corr_Nuc_2D_unpolarized[i] += v * SANSData.S_Nuc_2D_unpolarized[k] * c;
            SANSData.Corr_Mag_2D_unpolarized[i] += v * SANSData.S_Mag_2D_unpolarized[k] * c;
            SANSData.Corr_NucMag_2D[i] += v * SANSData.S_NucMag_2D[k] * c;
            SANSData.Corr_Mag_2D_polarized[i] += v * SANSData.S_Mag_2D_polarized[k] * c;
            SANSData.Corr_Mag_2D_chiral[i] += v * SANSData.S_Mag_2D_chiral[k] * c;
        }
    }
}





__global__
void Atomistic_MagSANS_Kernel_dilute(MagnetizationData MagData,\
                                     ScatteringData SANSData){

    // Input information:
    // N     : number of atoms
    // L     : number of points in Fourier space L = N_q*N_theta
    // K     : number of particles
    // x     : x-real-space position data in units of nano-meters
    // y     : y-real-space position data in units of nano-meters
    // z     : z-real-space position data in units of nano-meters
    // mx    : mx-real-space magnetic moment data in units of Bohr-Magneton
    // my    : my-real-space magnetic moment data in units of Bohr-Magneton
    // mz    : mz-real-space magnetci moment data in units of Bohr-Magneton
    // qy    : qy-Fourier-space coordinate in units of inverse nano-meters
    // qz    : qz-Fourier-space coordinate in units of inverse nano-meters
    // theta : theta-angle in Fourier space [theta = arctan2(qz, qy)] in radiant

    // output information:
    // Gxx_real: real-part of the xx-component of the Fourier-space correlation function of the magnetization
    // Gxx_imag: imaginary-part of the xx-component of the Fourier-space correlation function of the magnetization ...

    // dSigma_dOmega_M_unpolarized  : magnetic unpolarized SANS cross section
    // dSigma_dOmega_M_spin_flip    : magnetic spin-flip SANS cross section (sum of pm and mp )/2
    // dSigma_dOmega_M_chiral       : magnetic chiral SANS cross section (difference of pm and mp)/2
    // dSigma_dOmega_M_spin_flip_pm : magnetic pm spin-flip SANS cross section
    // dSigma_dOmega_M_spin_flip_mp : magnetic mp spin-flip SANS cross section

	unsigned long int L = (*SANSData.N_q) * (*SANSData.N_theta);
	unsigned long int N_avg = *MagData.N_avg;
	unsigned long int K = *MagData.K;
	unsigned long int W = *MagData.TotalAtomNumber;

	//float v = 1.0/((float)  (*MagData.K)) * pow(1.0/((float) (*MagData.N)), 2); // pre factor
	//float v = 1.0/((float) W);
	float v = 1.0/((float) W) * 1.0/((float) N_avg);
	int i = blockIdx.x * blockDim.x + threadIdx.x;

    float Px = SANSData.Polarization[0];
    float Py = SANSData.Polarization[1];
    float Pz = SANSData.Polarization[2];

	float mx_real = 0.0;
	float mx_imag = 0.0;
	float my_real = 0.0;
	float my_imag = 0.0;
	float mz_real = 0.0;
	float mz_imag = 0.0;

	float Mx_real = 0.0;
	float Mx_imag = 0.0;
	float My_real = 0.0;
	float My_imag = 0.0;
	float Mz_real = 0.0;
	float Mz_imag = 0.0;

	float Qx_real = 0.0;
	float Qx_imag = 0.0;
	float Qy_real = 0.0;
	float Qy_imag = 0.0;
	float Qz_real = 0.0;
	float Qz_imag = 0.0;


  //float X = 0.0;
	float Y = 0.0;
	float Z = 0.0;

	float Psi = 0.0;
	float cos_val = 0.0;
	float sin_val = 0.0;

	float cos_theta = 0.0;
	float sin_theta = 0.0;

	unsigned long int N_cum = 0;
	unsigned long int N_act = 0;
 
	if(i < L){
		for(int k=0; k < K; k++){

			mx_real = 0.0;
			mx_imag = 0.0;
			my_real = 0.0; 
			my_imag = 0.0;
			mz_real = 0.0; 
			mz_imag = 0.0;

        	Mx_real = 0.0;
        	Mx_imag = 0.0;
        	My_real = 0.0;
        	My_imag = 0.0;
        	Mz_real = 0.0;
        	Mz_imag = 0.0;

			Qx_real = 0.0;
			Qx_imag = 0.0;
			Qy_real = 0.0;
			Qy_imag = 0.0;
			Qz_real = 0.0;
			Qz_imag = 0.0;

			N_cum = MagData.N_cum[k];
			N_act = MagData.N_act[k];

        	for(int l=0; l < N_act; l++){
				// atomic position composition
				//X = MagData.RotMat[0] * (MagData.x[l+k*N] + StructData.x[k]) \
                //  + MagData.RotMat[3] * (MagData.y[l+k*N] + StructData.y[k]) \
				// + MagData.RotMat[6] * (MagData.z[l+k*N] + StructData.z[k]);
            	Y = MagData.RotMat[1] * MagData.x[l+N_cum] \
            	  + MagData.RotMat[4] * MagData.y[l+N_cum] \
            	  + MagData.RotMat[7] * MagData.z[l+N_cum];
            	Z = MagData.RotMat[2] * MagData.x[l+N_cum] \
            	  + MagData.RotMat[5] * MagData.y[l+N_cum] \
            	  + MagData.RotMat[8] * MagData.z[l+N_cum];

				// phase function
				Psi = Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i];

				// cosine and sine values
				cos_val = cos(Psi);
				sin_val = sin(Psi);

				// cosine and sine summations
            	mx_real += MagData.mx[l+N_cum] * cos_val;
            	mx_imag -= MagData.mx[l+N_cum] * sin_val;
            	my_real += MagData.my[l+N_cum] * cos_val;
            	my_imag -= MagData.my[l+N_cum] * sin_val;
            	mz_real += MagData.mz[l+N_cum] * cos_val;
            	mz_imag -= MagData.mz[l+N_cum] * sin_val;

			}

			Mx_real = MagData.RotMat[0] * mx_real + MagData.RotMat[3] * my_real + MagData.RotMat[6] * mz_real;
			My_real = MagData.RotMat[1] * mx_real + MagData.RotMat[4] * my_real + MagData.RotMat[7] * mz_real;
			Mz_real = MagData.RotMat[2] * mx_real + MagData.RotMat[5] * my_real + MagData.RotMat[8] * mz_real;

			Mx_imag = MagData.RotMat[0] * mx_imag + MagData.RotMat[3] * my_imag + MagData.RotMat[6] * mz_imag;
			My_imag = MagData.RotMat[1] * mx_imag + MagData.RotMat[4] * my_imag + MagData.RotMat[7] * mz_imag;
			Mz_imag = MagData.RotMat[2] * mx_imag + MagData.RotMat[5] * my_imag + MagData.RotMat[8] * mz_imag;


			cos_theta = cos(SANSData.theta_2D[i]);
			sin_theta = sin(SANSData.theta_2D[i]);

			// real-parts of the Halpern-Johnson vector
			Qx_real = (-Mx_real);
			Qy_real = (Mz_real * sin_theta - My_real * cos_theta) * cos_theta;
			Qz_real = (My_real * cos_theta - Mz_real * sin_theta) * sin_theta;

			// imaginary-parts of the Halpern-Johnson vector
			Qx_imag = (-Mx_imag);
			Qy_imag = (Mz_imag * sin_theta - My_imag * cos_theta) * cos_theta;
			Qz_imag = (My_imag * cos_theta - Mz_imag * sin_theta) * sin_theta;


			// nuclear SANS cross section projected in (qz, qy)-plane
			SANSData.S_Nuc_2D_unpolarized[i] += 0.0;

			// unpolarized magnetic SANS cross section projected in (qz, qy)-plane
			SANSData.S_Mag_2D_unpolarized[i] += v * (Qx_real * Qx_real + Qx_imag * Qx_imag) \
											  + v * (Qy_real * Qy_real + Qy_imag * Qy_imag) \
											  + v * (Qz_real * Qz_real + Qz_imag * Qz_imag);

			// nuclear magnetic interference SANS cross section projected in (qz, qy)-plane
			SANSData.S_NucMag_2D[i] += 0.0;

			// polarized magnetic SANS cross section projected in the (qz, qy)-plane
			SANSData.S_Mag_2D_polarized[i] += v * pow(Px, 2) * (Qx_real * Qx_real + Qx_imag * Qx_imag) \
										    + v * pow(Py, 2) * (Qy_real * Qy_real + Qy_imag * Qy_imag) \
										    + v * pow(Pz, 2) * (Qz_real * Qz_real + Qz_imag * Qz_imag) \
										    + v * 2.0 * Px * Py * (Qx_real * Qy_real + Qx_imag * Qy_imag) \
										    + v * 2.0 * Px * Pz * (Qx_real * Qz_real + Qx_imag * Qz_imag) \
										    + v * 2.0 * Py * Pz * (Qy_real * Qz_real + Qy_imag * Qz_imag);

			// chiral magnetic SANS cross section in (qz, qy)-plane
			SANSData.S_Mag_2D_chiral[i] += v * 2.0 * Px * (Qy_imag * Qz_real - Qz_imag * Qy_real) \
                                         + v * 2.0 * Py * (Qz_imag * Qx_real - Qx_imag * Qz_real) \
                                         + v * 2.0 * Pz * (Qx_imag * Qy_real - Qy_imag * Qx_real);


			SANSData.Gxx_real[i] += v*(Mx_real * Mx_real + Mx_imag * Mx_imag);
			SANSData.Gxx_imag[i] += 0.0;

			SANSData.Gyy_real[i] += v*(My_real * My_real + My_imag * My_imag);
			SANSData.Gyy_imag[i] += 0.0;

			SANSData.Gzz_real[i] += v*(Mz_real * Mz_real + Mz_imag * Mz_imag);
			SANSData.Gzz_imag[i] += 0.0;

			SANSData.Gxy_real[i] += v*(Mx_real * My_real + Mx_imag * My_imag);
			SANSData.Gxy_imag[i] += v*(Mx_imag * My_real - Mx_real * My_imag);

			SANSData.Gyx_real[i] =  SANSData.Gxy_real[i];
			SANSData.Gyx_imag[i] = -SANSData.Gxy_imag[i];

			SANSData.Gxz_real[i] += v*(Mx_real * Mz_real + Mx_imag * Mz_imag);
			SANSData.Gxz_imag[i] += v*(Mx_imag * Mz_real - Mx_real * Mz_imag);

			SANSData.Gzx_real[i] =  SANSData.Gxz_real[i];
			SANSData.Gzx_imag[i] = -SANSData.Gxy_imag[i];

			SANSData.Gyz_real[i] += v*(My_real * Mz_real + My_imag * Mz_imag);
			SANSData.Gyz_imag[i] += v*(My_imag * Mz_real - My_real * Mz_imag);

			SANSData.Gzx_real[i] =  SANSData.Gyz_real[i];
			SANSData.Gzx_imag[i] = -SANSData.Gyz_imag[i];

	}

	}
}



__global__
void Atomistic_NucSANS_Kernel_dilute(NuclearData NucData,\
                                     ScatteringData SANSData){

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

	//float v = 1.0/((float)  (*NucData.K)) * pow(1.0/((float) (*NucData.N)), 2); // pre factor
	float v = 1.0/((float) W) * 1.0/((float) N_avg);;

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
 
	if(i < L){
		for(int k=0; k< (*NucData.K); k++){

        	Nuc_real = 0.0;
        	Nuc_imag = 0.0;

			N_cum = NucData.N_cum[k];
			N_act = NucData.N_act[k];

        	for(int l=0; l < N_act; l++){
				// atomic position composition
				//X = MagData.RotMat[0] * MagData.x[l+k*N] \
                //  + MagData.RotMat[3] * MagData.y[l+k*N] \
				// + MagData.RotMat[6] * MagData.z[l+k*N];
            	Y = NucData.RotMat[1] * NucData.x[l+N_cum] \
            	  + NucData.RotMat[4] * NucData.y[l+N_cum] \
            	  + NucData.RotMat[7] * NucData.z[l+N_cum];
            	Z = NucData.RotMat[2] * NucData.x[l+N_cum] \
            	  + NucData.RotMat[5] * NucData.y[l+N_cum] \
            	  + NucData.RotMat[8] * NucData.z[l+N_cum];

				// phase function
				Psi = Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i];

				// cosine and sine values
				cos_val = cos(Psi);
				sin_val = sin(Psi);

				// cosine- and sine-summations
				Nuc_real += NucData.Nuc[l+N_cum] * cos_val;
				Nuc_imag -= NucData.Nuc[l+N_cum] * sin_val;
            	
			}

			// nuclear SANS cross section projected in (qz, qy)-plane
			SANSData.S_Nuc_2D_unpolarized[i] += v * (Nuc_real * Nuc_real + Nuc_imag * Nuc_imag);
			
		}

	}
}



__global__
void Atomistic_NuMagSANS_Kernel_dilute(NuclearData NucData,\
									   MagnetizationData MagData,\
                                       ScatteringData SANSData){

    // Input information:
    // N     : number of atoms
    // L     : number of points in Fourier space L = N_q*N_theta
    // K     : number of particles
   // x     : x-real-space position data in units of nano-meters
    // y     : y-real-space position data in units of nano-meters
   // z     : z-real-space position data in units of nano-meters
   // mx    : mx-real-space magnetic moment data in units of Bohr-Magneton
   // my    : my-real-space magnetic moment data in units of Bohr-Magneton
    // mz    : mz-real-space magnetci moment data in units of Bohr-Magneton
   // qy    : qy-Fourier-space coordinate in units of inverse nano-meters
   // qz    : qz-Fourier-space coordinate in units of inverse nano-meters
   // theta : theta-angle in Fourier space [theta = arctan2(qz, qy)] in radiant

   // output information:
   // Gxx_real: real-part of the xx-component of the Fourier-space correlation function of the magnetization
   // Gxx_imag: imaginary-part of the xx-component of the Fourier-space correlation function of the magnetization ...

   // dSigma_dOmega_M_unpolarized  : magnetic unpolarized SANS cross section
   // dSigma_dOmega_M_spin_flip    : magnetic spin-flip SANS cross section (sum of pm and mp )/2
   // dSigma_dOmega_M_chiral       : magnetic chiral SANS cross section (difference of pm and mp)/2
   // dSigma_dOmega_M_spin_flip_pm : magnetic pm spin-flip SANS cross section
   // dSigma_dOmega_M_spin_flip_mp : magnetic mp spin-flip SANS cross section

	unsigned long int L = (*SANSData.N_q) * (*SANSData.N_theta);
	unsigned long int N_avg = *MagData.N_avg;
	unsigned long int K = *MagData.K;
	unsigned long int W = *MagData.TotalAtomNumber;

	//float v = (1.0/((float) K)) * pow(1.0/((float) N), 2); // pre factor
	float v = 1.0/((float) W) * 1.0/((float) N_avg);
	int i = blockIdx.x * blockDim.x + threadIdx.x;

    float Px = SANSData.Polarization[0];
    float Py = SANSData.Polarization[1];
    float Pz = SANSData.Polarization[2];

	float mx_real = 0.0;
	float mx_imag = 0.0;
	float my_real = 0.0;
	float my_imag = 0.0;
	float mz_real = 0.0;
	float mz_imag = 0.0;

	float Mx_real = 0.0;
	float Mx_imag = 0.0;
	float My_real = 0.0;
	float My_imag = 0.0;
	float Mz_real = 0.0;
	float Mz_imag = 0.0;

	float Qx_real = 0.0;
	float Qx_imag = 0.0;
	float Qy_real = 0.0;
	float Qy_imag = 0.0;
	float Qz_real = 0.0;
	float Qz_imag = 0.0;

	float nuc_real = 0.0;
	float nuc_imag = 0.0;
   // float X = 0.0;
	float Y = 0.0;
	float Z = 0.0;

	float Psi = 0.0;

	float cos_val = 0.0;
	float sin_val = 0.0;

	float cos_theta = 0.0;
	float sin_theta = 0.0;

	unsigned long int N_cum = 0;
	unsigned long int N_act = 0;
 
	if(i < L){
		for(int k=0; k < K; k++){

			mx_real = 0.0;
			mx_imag = 0.0;
			my_real = 0.0;
			my_imag = 0.0;
			mz_real = 0.0;
			mz_imag = 0.0;
			
        	Mx_real = 0.0;
        	Mx_imag = 0.0;
        	My_real = 0.0;
        	My_imag = 0.0;
        	Mz_real = 0.0;
        	Mz_imag = 0.0;

			Qx_real = 0.0;
			Qx_imag = 0.0;
			Qy_real = 0.0;
			Qy_imag = 0.0;
			Qz_real = 0.0;
			Qz_imag = 0.0;

        	nuc_real = 0.0;
        	nuc_imag = 0.0;

			N_cum = MagData.N_cum[k];
			N_act = MagData.N_act[k];

        	for(int l=0; l < N_act; l++){
				// atomic position composition
				//X = MagData.RotMat[0] * MagData.x[l+k*N] \
                //  + MagData.RotMat[3] * MagData.y[l+k*N] \
                //  + MagData.RotMat[6] * MagData.z[l+k*N];
            	Y = MagData.RotMat[1] * MagData.x[l+N_cum] \
            	  + MagData.RotMat[4] * MagData.y[l+N_cum] \
            	  + MagData.RotMat[7] * MagData.z[l+N_cum];
            	Z = MagData.RotMat[2] * MagData.x[l+N_cum] \
            	  + MagData.RotMat[5] * MagData.y[l+N_cum] \
            	  + MagData.RotMat[8] * MagData.z[l+N_cum];

				// phase function
				Psi = Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i];

				// cosine and sine values
				cos_val = cos(Psi);
				sin_val = sin(Psi);

				// cosine and sine summations
				nuc_real += NucData.Nuc[l+N_cum] * cos_val;
				nuc_imag -= NucData.Nuc[l+N_cum] * sin_val;

            	mx_real += MagData.mx[l+N_cum] * cos_val;
            	mx_imag -= MagData.mx[l+N_cum] * sin_val;
            	my_real += MagData.my[l+N_cum] * cos_val;
            	my_imag -= MagData.my[l+N_cum] * sin_val;
            	mz_real += MagData.mz[l+N_cum] * cos_val;
            	mz_imag -= MagData.mz[l+N_cum] * sin_val;

			}

			// rotations of the magnetization fourier components real part
			Mx_real = MagData.RotMat[0] * mx_real + MagData.RotMat[3] * my_real + MagData.RotMat[6] * mz_real;
			My_real = MagData.RotMat[1] * mx_real + MagData.RotMat[4] * my_real + MagData.RotMat[7] * mz_real;
			Mz_real = MagData.RotMat[2] * mx_real + MagData.RotMat[5] * my_real + MagData.RotMat[8] * mz_real;

			// rotations of the magnetization fourier components imaginary part
			Mx_imag = MagData.RotMat[0] * mx_imag + MagData.RotMat[3] * my_imag + MagData.RotMat[6] * mz_imag;
			My_imag = MagData.RotMat[1] * mx_imag + MagData.RotMat[4] * my_imag + MagData.RotMat[7] * mz_imag;
			Mz_imag = MagData.RotMat[2] * mx_imag + MagData.RotMat[5] * my_imag + MagData.RotMat[8] * mz_imag;


			cos_theta = cos(SANSData.theta_2D[i]);
			sin_theta = sin(SANSData.theta_2D[i]);

			// real-parts of the Halpern-Johnson vector
			Qx_real = (-Mx_real);
			Qy_real = (Mz_real * sin_theta - My_real * cos_theta) * cos_theta;
			Qz_real = (My_real * cos_theta - Mz_real * sin_theta) * sin_theta;

			// imaginary-parts of the Halpern-Johnson vector
			Qx_imag = (-Mx_imag);
			Qy_imag = (Mz_imag * sin_theta - My_imag * cos_theta) * cos_theta;
			Qz_imag = (My_imag * cos_theta - Mz_imag * sin_theta) * sin_theta;

			// nuclear SANS cross section projected in (qz, qy)-plane
			SANSData.S_Nuc_2D_unpolarized[i] += v * (nuc_real * nuc_real + nuc_imag * nuc_imag);

			// unpolarized magnetic SANS cross section projected in (qz, qy)-plane
			SANSData.S_Mag_2D_unpolarized[i] += v * (Qx_real * Qx_real + Qx_imag * Qx_imag) \
											  + v * (Qy_real * Qy_real + Qy_imag * Qy_imag) \
											  + v * (Qz_real * Qz_real + Qz_imag * Qz_imag);

			// nuclear magnetic interference SANS cross section projected in (qz, qy)-plane
			SANSData.S_NucMag_2D[i] += 2.0 * v * Px * (nuc_real * Qx_real + nuc_imag * Qx_imag) \
									 + 2.0 * v * Py * (nuc_real * Qy_real + nuc_imag * Qy_imag) \
									 + 2.0 * v * Pz * (nuc_real * Qz_real + nuc_imag * Qz_imag);

			// polarized magnetic SANS cross section projected in the (qz, qy)-plane
			SANSData.S_Mag_2D_polarized[i] += v * pow(Px, 2) * (Qx_real * Qx_real + Qx_imag * Qx_imag) \
										    + v * pow(Py, 2) * (Qy_real * Qy_real + Qy_imag * Qy_imag) \
										    + v * pow(Pz, 2) * (Qz_real * Qz_real + Qz_imag * Qz_imag) \
										    + v * 2.0 * Px * Py * (Qx_real * Qy_real + Qx_imag * Qy_imag) \
										    + v * 2.0 * Px * Pz * (Qx_real * Qz_real + Qx_imag * Qz_imag) \
										    + v * 2.0 * Py * Pz * (Qy_real * Qz_real + Qy_imag * Qz_imag);

			// chiral magnetic SANS cross section in (qz, qy)-plane
			SANSData.S_Mag_2D_chiral[i] += v * 2.0 * Px * (Qy_imag * Qz_real - Qz_imag * Qy_real) \
                                         + v * 2.0 * Py * (Qz_imag * Qx_real - Qx_imag * Qz_real) \
                                         + v * 2.0 * Pz * (Qx_imag * Qy_real - Qy_imag * Qx_real);


			SANSData.Gxx_real[i] += v*(Mx_real * Mx_real + Mx_imag * Mx_imag);
			SANSData.Gxx_imag[i] += 0.0;

			SANSData.Gyy_real[i] += v*(My_real * My_real + My_imag * My_imag);
			SANSData.Gyy_imag[i] += 0.0;

			SANSData.Gzz_real[i] += v*(Mz_real * Mz_real + Mz_imag * Mz_imag);
			SANSData.Gzz_imag[i] += 0.0;

			SANSData.Gxy_real[i] += v*(Mx_real * My_real + Mx_imag * My_imag);
			SANSData.Gxy_imag[i] += v*(Mx_imag * My_real - Mx_real * My_imag);

			SANSData.Gyx_real[i] =  SANSData.Gxy_real[i];
			SANSData.Gyx_imag[i] = -SANSData.Gxy_imag[i];

			SANSData.Gxz_real[i] += v*(Mx_real * Mz_real + Mx_imag * Mz_imag);
			SANSData.Gxz_imag[i] += v*(Mx_imag * Mz_real - Mx_real * Mz_imag);

			SANSData.Gzx_real[i] =  SANSData.Gxz_real[i];
			SANSData.Gzx_imag[i] = -SANSData.Gxy_imag[i];

			SANSData.Gyz_real[i] += v*(My_real * Mz_real + My_imag * Mz_imag);
			SANSData.Gyz_imag[i] += v*(My_imag * Mz_real - My_real * Mz_imag);

			SANSData.Gzx_real[i] =  SANSData.Gyz_real[i];
			SANSData.Gzx_imag[i] = -SANSData.Gyz_imag[i];

	}

	}
}






__global__
void Atomistic_MagSANS_Kernel(MagnetizationData MagData,\
							  StructureData StructData, \
							  ScatteringData SANSData){

    // Input information:
    // N     : number of atoms
    // L     : number of points in Fourier space L = N_q*N_theta
    // K     : number of particles
   // x     : x-real-space position data in units of nano-meters
    // y     : y-real-space position data in units of nano-meters
   // z     : z-real-space position data in units of nano-meters
   // mx    : mx-real-space magnetic moment data in units of Bohr-Magneton
   // my    : my-real-space magnetic moment data in units of Bohr-Magneton
    // mz    : mz-real-space magnetci moment data in units of Bohr-Magneton
   // qy    : qy-Fourier-space coordinate in units of inverse nano-meters
   // qz    : qz-Fourier-space coordinate in units of inverse nano-meters
   // theta : theta-angle in Fourier space [theta = arctan2(qz, qy)] in radiant

   // output information:
   // Gxx_real: real-part of the xx-component of the Fourier-space correlation function of the magnetization
   // Gxx_imag: imaginary-part of the xx-component of the Fourier-space correlation function of the magnetization ...

   // dSigma_dOmega_M_unpolarized  : magnetic unpolarized SANS cross section
   // dSigma_dOmega_M_spin_flip    : magnetic spin-flip SANS cross section (sum of pm and mp )/2
   // dSigma_dOmega_M_chiral       : magnetic chiral SANS cross section (difference of pm and mp)/2
   // dSigma_dOmega_M_spin_flip_pm : magnetic pm spin-flip SANS cross section
   // dSigma_dOmega_M_spin_flip_mp : magnetic mp spin-flip SANS cross section

	unsigned long int L = (*SANSData.N_q) * (*SANSData.N_theta);
	unsigned long int N_avg = *MagData.N_avg;
	unsigned long int W = *MagData.TotalAtomNumber;

	//float v = 1.0/((float)  (*MagData.K)) * pow(1.0/((float) (*MagData.N)), 2); // pre factor
	float v =  1.0/((float) W) * 1.0/((float) N_avg);

	int i = blockIdx.x * blockDim.x + threadIdx.x;

    float Px = SANSData.Polarization[0];
    float Py = SANSData.Polarization[1];
    float Pz = SANSData.Polarization[2];

	float mx_real = 0.0;
	float mx_imag = 0.0;
	float my_real = 0.0;
	float my_imag = 0.0;
	float mz_real = 0.0;
	float mz_imag = 0.0;

	float Mx_real = 0.0;
	float Mx_imag = 0.0;
	float My_real = 0.0;
	float My_imag = 0.0;
	float Mz_real = 0.0;
	float Mz_imag = 0.0;

	float Qx_real = 0.0;
	float Qx_imag = 0.0;
	float Qy_real = 0.0;
	float Qy_imag = 0.0;
	float Qz_real = 0.0;
	float Qz_imag = 0.0;

  //float X = 0.0;
	float Y = 0.0;
	float Z = 0.0;

	float Psi = 0.0;
	float cos_val = 0.0;
	float sin_val = 0.0;

	float cos_theta = 0.0;
	float sin_theta = 0.0;

	unsigned long int N_cum = 0;
	unsigned long int N_act = 0;

	if(i < L){
		for(int k=0; k< (*MagData.K); k++){

			mx_real = 0.0;
			mx_imag = 0.0;
			my_real = 0.0;
			my_imag = 0.0;
			mz_real = 0.0;
			mz_imag = 0.0;

			N_cum = MagData.N_cum[k];
			N_act = MagData.N_act[k];

        	for(int l=0; l < N_act; l++){
				// atomic position composition
				//X = MagData.RotMat[0] * (MagData.x[l+k*N] + StructData.x[k]) \
                //  + MagData.RotMat[3] * (MagData.y[l+k*N] + StructData.y[k]) \
				// + MagData.RotMat[6] * (MagData.z[l+k*N] + StructData.z[k]);
            	Y = MagData.RotMat[1] * (MagData.x[l+N_cum] + StructData.x[k]) \
            	  + MagData.RotMat[4] * (MagData.y[l+N_cum] + StructData.y[k]) \
            	  + MagData.RotMat[7] * (MagData.z[l+N_cum] + StructData.z[k]);
            	Z = MagData.RotMat[2] * (MagData.x[l+N_cum] + StructData.x[k]) \
            	  + MagData.RotMat[5] * (MagData.y[l+N_cum] + StructData.y[k]) \
            	  + MagData.RotMat[8] * (MagData.z[l+N_cum] + StructData.z[k]);

				// phase function
				Psi = Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i];

				// cosine and sine values
				cos_val = cos(Psi);
				sin_val = sin(Psi);

				// cosine and sine summations
            	mx_real += MagData.mx[l+N_cum] * cos_val;
            	mx_imag -= MagData.mx[l+N_cum] * sin_val;
            	my_real += MagData.my[l+N_cum] * cos_val;
            	my_imag -= MagData.my[l+N_cum] * sin_val;
            	mz_real += MagData.mz[l+N_cum] * cos_val;
            	mz_imag -= MagData.mz[l+N_cum] * sin_val;

			}

			Mx_real += MagData.RotMat[0] * mx_real + MagData.RotMat[3] * my_real + MagData.RotMat[6] * mz_real;
			My_real += MagData.RotMat[1] * mx_real + MagData.RotMat[4] * my_real + MagData.RotMat[7] * mz_real;
			Mz_real += MagData.RotMat[2] * mx_real + MagData.RotMat[5] * my_real + MagData.RotMat[8] * mz_real;

			Mx_imag += MagData.RotMat[0] * mx_imag + MagData.RotMat[3] * my_imag + MagData.RotMat[6] * mz_imag;
			My_imag += MagData.RotMat[1] * mx_imag + MagData.RotMat[4] * my_imag + MagData.RotMat[7] * mz_imag;
			Mz_imag += MagData.RotMat[2] * mx_imag + MagData.RotMat[5] * my_imag + MagData.RotMat[8] * mz_imag;

		}

		cos_theta = cos(SANSData.theta_2D[i]);
		sin_theta = sin(SANSData.theta_2D[i]);

		// real-parts of the Halpern-Johnson vector
		Qx_real = (-Mx_real);
		Qy_real = (Mz_real * sin_theta - My_real * cos_theta) * cos_theta;
		Qz_real = (My_real * cos_theta - Mz_real * sin_theta) * sin_theta;

		// imaginary-parts of the Halpern-Johnson vector
		Qx_imag = (-Mx_imag);
		Qy_imag = (Mz_imag * sin_theta - My_imag * cos_theta) * cos_theta;
		Qz_imag = (My_imag * cos_theta - Mz_imag * sin_theta) * sin_theta;


		// nuclear SANS cross section projected in (qz, qy)-plane
		SANSData.S_Nuc_2D_unpolarized[i] = 0.0;

		// unpolarized magnetic SANS cross section projected in (qz, qy)-plane
		SANSData.S_Mag_2D_unpolarized[i] = v * (Qx_real * Qx_real + Qx_imag * Qx_imag) \
										 + v * (Qy_real * Qy_real + Qy_imag * Qy_imag) \
										 + v * (Qz_real * Qz_real + Qz_imag * Qz_imag);

		// nuclear magnetic interference SANS cross section projected in (qz, qy)-plane
		SANSData.S_NucMag_2D[i] = 0.0;

		// polarized magnetic SANS cross section projected in the (qz, qy)-plane
		SANSData.S_Mag_2D_polarized[i] = v * pow(Px, 2) * (Qx_real * Qx_real + Qx_imag * Qx_imag) \
									   + v * pow(Py, 2) * (Qy_real * Qy_real + Qy_imag * Qy_imag) \
									   + v * pow(Pz, 2) * (Qz_real * Qz_real + Qz_imag * Qz_imag) \
									   + v * 2.0 * Px * Py * (Qx_real * Qy_real + Qx_imag * Qy_imag) \
									   + v * 2.0 * Px * Pz * (Qx_real * Qz_real + Qx_imag * Qz_imag) \
									   + v * 2.0 * Py * Pz * (Qy_real * Qz_real + Qy_imag * Qz_imag);

		// chiral magnetic SANS cross section in (qz, qy)-plane
		SANSData.S_Mag_2D_chiral[i] = v * 2.0 * Px * (Qy_imag * Qz_real - Qz_imag * Qy_real) \
									+ v * 2.0 * Py * (Qz_imag * Qx_real - Qx_imag * Qz_real) \
									+ v * 2.0 * Pz * (Qx_imag * Qy_real - Qy_imag * Qx_real);


		SANSData.Gxx_real[i] = v*(Mx_real * Mx_real + Mx_imag * Mx_imag);
		SANSData.Gxx_imag[i] = 0.0;

		SANSData.Gyy_real[i] = v*(My_real * My_real + My_imag * My_imag);
		SANSData.Gyy_imag[i] = 0.0;

		SANSData.Gzz_real[i] = v*(Mz_real * Mz_real + Mz_imag * Mz_imag);
		SANSData.Gzz_imag[i] = 0.0;

		SANSData.Gxy_real[i] = v*(Mx_real * My_real + Mx_imag * My_imag);
		SANSData.Gxy_imag[i] = v*(Mx_imag * My_real - Mx_real * My_imag);

		SANSData.Gyx_real[i] =  SANSData.Gxy_real[i];
		SANSData.Gyx_imag[i] = -SANSData.Gxy_imag[i];

		SANSData.Gxz_real[i] = v*(Mx_real * Mz_real + Mx_imag * Mz_imag);
		SANSData.Gxz_imag[i] = v*(Mx_imag * Mz_real - Mx_real * Mz_imag);

		SANSData.Gzx_real[i] =  SANSData.Gxz_real[i];
		SANSData.Gzx_imag[i] = -SANSData.Gxy_imag[i];

		SANSData.Gyz_real[i] = v*(My_real * Mz_real + My_imag * Mz_imag);
		SANSData.Gyz_imag[i] = v*(My_imag * Mz_real - My_real * Mz_imag);

		SANSData.Gzx_real[i] =  SANSData.Gyz_real[i];
		SANSData.Gzx_imag[i] = -SANSData.Gyz_imag[i];

	}
}





__global__
void Atomistic_NucSANS_Kernel(NuclearData NucData,\
							  StructureData StructData, \
							  ScatteringData SANSData){

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
	unsigned long int K = *NucData.K;
	unsigned long int W = *NucData.TotalAtomNumber;

	//float v = 1.0/((float)  (*NucData.K)) * pow(1.0/((float) (*NucData.N)), 2); // pre factor
	float v = 1.0/((float) W) * 1.0/((float) N_avg);

	int i = blockIdx.x * blockDim.x + threadIdx.x;

	float nuc_real = 0.0;
	float nuc_imag = 0.0;

	float Nuc_real = 0.0;
	float Nuc_imag = 0.0;

	float Y = 0.0;
	float Z = 0.0;

	float Psi = 0.0;

	float cos_val = 0.0;
	float sin_val = 0.0;

	unsigned long int N_cum = 0;
	unsigned long int N_act = 0;

	if(i < L){
		for(int k=0; k < K; k++){

        	nuc_real = 0.0;
        	nuc_imag = 0.0;

			N_cum = NucData.N_cum[k];
			N_act = NucData.N_act[k];

        	for(int l=0; l < N_act; l++){
				// atomic position composition
				//X = MagData.RotMat[0] * (MagData.x[l+k*N] + StructData.x[k]) \
                //  + MagData.RotMat[3] * (MagData.y[l+k*N] + StructData.y[k]) \
				// + MagData.RotMat[6] * (MagData.z[l+k*N] + StructData.z[k]);
            	Y = NucData.RotMat[1] * (NucData.x[l+N_cum] + StructData.x[k]) \
            	  + NucData.RotMat[4] * (NucData.y[l+N_cum] + StructData.y[k]) \
            	  + NucData.RotMat[7] * (NucData.z[l+N_cum] + StructData.z[k]);
            	Z = NucData.RotMat[2] * (NucData.x[l+N_cum] + StructData.x[k]) \
            	  + NucData.RotMat[5] * (NucData.y[l+N_cum] + StructData.y[k]) \
            	  + NucData.RotMat[8] * (NucData.z[l+N_cum] + StructData.z[k]);

				// phase function
				Psi = Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i];

				// cosine and sine values
				cos_val = cos(Psi);
				sin_val = sin(Psi);

				nuc_real += NucData.Nuc[l+N_cum] * cos_val;
				nuc_imag -= NucData.Nuc[l+N_cum] * sin_val;

			}

			Nuc_real += nuc_real;
			Nuc_imag += nuc_imag;

		}

		// nuclear SANS cross section projected in (qz, qy)-plane
		SANSData.S_Nuc_2D_unpolarized[i] = v * (Nuc_real * Nuc_real + Nuc_imag * Nuc_imag);

	}
}





__global__
void Atomistic_NuMagSANS_Kernel(NuclearData NucData, \
							    MagnetizationData MagData, \
							    StructureData StructData, \
								ScatteringData SANSData){

    // Input information:
    // N     : number of atoms
    // L     : number of points in Fourier space L = N_q*N_theta
    // K     : number of particles
    // x     : x-real-space position data in units of nano-meters
    // y     : y-real-space position data in units of nano-meters
    // z     : z-real-space position data in units of nano-meters
    // mx    : mx-real-space magnetic moment data in units of Bohr-Magneton
    // my    : my-real-space magnetic moment data in units of Bohr-Magneton
    // mz    : mz-real-space magnetci moment data in units of Bohr-Magneton
    // qy    : qy-Fourier-space coordinate in units of inverse nano-meters
    // qz    : qz-Fourier-space coordinate in units of inverse nano-meters
    // theta : theta-angle in Fourier space [theta = arctan2(qz, qy)] in radiant

    // output information:
    // Gxx_real: real-part of the xx-component of the Fourier-space correlation function of the magnetization
    // Gxx_imag: imaginary-part of the xx-component of the Fourier-space correlation function of the magnetization ...

    // dSigma_dOmega_M_unpolarized  : magnetic unpolarized SANS cross section
    // dSigma_dOmega_M_spin_flip    : magnetic spin-flip SANS cross section (sum of pm and mp )/2
    // dSigma_dOmega_M_chiral       : magnetic chiral SANS cross section (difference of pm and mp)/2
    // dSigma_dOmega_M_spin_flip_pm : magnetic pm spin-flip SANS cross section
    // dSigma_dOmega_M_spin_flip_mp : magnetic mp spin-flip SANS cross section

	unsigned long int L = (*SANSData.N_q) * (*SANSData.N_theta);
	unsigned long int N_avg = *MagData.N_avg;
	unsigned long int K = *MagData.K;
	unsigned long int W = *MagData.TotalAtomNumber;

	//float v = (1.0/((float) K)) * pow(1.0/((float) N), 2); // pre factor
	float v = 1.0/((float) W) * 1.0/((float) N_avg);

	int i = blockIdx.x * blockDim.x + threadIdx.x;

    float Px = SANSData.Polarization[0];
    float Py = SANSData.Polarization[1];
    float Pz = SANSData.Polarization[2];

	float mx_real = 0.0;
	float mx_imag = 0.0;
	float my_real = 0.0;
	float my_imag = 0.0;
	float mz_real = 0.0;
	float mz_imag = 0.0;

	float Mx_real = 0.0;
	float Mx_imag = 0.0;
	float My_real = 0.0;
	float My_imag = 0.0;
	float Mz_real = 0.0;
	float Mz_imag = 0.0;

	float Qx_real = 0.0;
	float Qx_imag = 0.0;
	float Qy_real = 0.0;
	float Qy_imag = 0.0;
	float Qz_real = 0.0;
	float Qz_imag = 0.0;

	float nuc_real = 0.0;
	float nuc_imag = 0.0;

	float Nuc_real = 0.0;
	float Nuc_imag = 0.0;
   // float X = 0.0;
	float Y = 0.0;
	float Z = 0.0;

	float Psi = 0.0;

	float cos_theta = 0.0;
	float sin_theta = 0.0;

	float cos_val = 0.0;
	float sin_val = 0.0;

	unsigned long int N_cum = 0;
	unsigned long int N_act = 0;


	if(i < L){
		for(int k=0; k < K; k++){

			mx_real = 0.0;
			mx_imag = 0.0;
			my_real = 0.0;
			my_imag = 0.0;
			mz_real = 0.0;
			mz_imag = 0.0;

        	nuc_real = 0.0;
        	nuc_imag = 0.0;

			N_cum = MagData.N_cum[k];
			N_act = MagData.N_act[k];

        	for(int l=0; l < N_act; l++){

				// atomic position composition
				//X = MagData.RotMat[0] * (MagData.x[l+k*N] + StructData.x[k]) \
                //  + MagData.RotMat[3] * (MagData.y[l+k*N] + StructData.y[k]) \
				// + MagData.RotMat[6] * (MagData.z[l+k*N] + StructData.z[k]);
            	Y = MagData.RotMat[1] * (MagData.x[l+N_cum] + StructData.x[k]) \
            	  + MagData.RotMat[4] * (MagData.y[l+N_cum] + StructData.y[k]) \
            	  + MagData.RotMat[7] * (MagData.z[l+N_cum] + StructData.z[k]);
            	Z = MagData.RotMat[2] * (MagData.x[l+N_cum] + StructData.x[k]) \
            	  + MagData.RotMat[5] * (MagData.y[l+N_cum] + StructData.y[k]) \
            	  + MagData.RotMat[8] * (MagData.z[l+N_cum] + StructData.z[k]);

				// phase function
				Psi = Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i];

				// cosine and sine values
				cos_val = cos(Psi);
				sin_val = sin(Psi);

				// cosine and sine summations
				nuc_real += NucData.Nuc[l+N_cum] * cos_val;
				nuc_imag -= NucData.Nuc[l+N_cum] * sin_val;

            	mx_real += MagData.mx[l+N_cum] * cos_val;
            	mx_imag -= MagData.mx[l+N_cum] * sin_val;
            	my_real += MagData.my[l+N_cum] * cos_val;
            	my_imag -= MagData.my[l+N_cum] * sin_val;
            	mz_real += MagData.mz[l+N_cum] * cos_val;
            	mz_imag -= MagData.mz[l+N_cum] * sin_val;
			}

			Nuc_real += nuc_real;
			Nuc_imag += nuc_imag;

			// rotations of the magnetization fourier components real part
			Mx_real += MagData.RotMat[0] * mx_real + MagData.RotMat[3] * my_real + MagData.RotMat[6] * mz_real;
			My_real += MagData.RotMat[1] * mx_real + MagData.RotMat[4] * my_real + MagData.RotMat[7] * mz_real;
			Mz_real += MagData.RotMat[2] * mx_real + MagData.RotMat[5] * my_real + MagData.RotMat[8] * mz_real;

			// rotations of the magnetization fourier components imaginary part
			Mx_imag += MagData.RotMat[0] * mx_imag + MagData.RotMat[3] * my_imag + MagData.RotMat[6] * mz_imag;
			My_imag += MagData.RotMat[1] * mx_imag + MagData.RotMat[4] * my_imag + MagData.RotMat[7] * mz_imag;
			Mz_imag += MagData.RotMat[2] * mx_imag + MagData.RotMat[5] * my_imag + MagData.RotMat[8] * mz_imag;

		}

		cos_theta = cos(SANSData.theta_2D[i]);
		sin_theta = sin(SANSData.theta_2D[i]);

		// real-parts of the Halpern-Johnson vector
		Qx_real = (-Mx_real);
		Qy_real = (Mz_real * sin_theta - My_real * cos_theta) * cos_theta;
		Qz_real = (My_real * cos_theta - Mz_real * sin_theta) * sin_theta;

		// imaginary-parts of the Halpern-Johnson vector
		Qx_imag = (-Mx_imag);
		Qy_imag = (Mz_imag * sin_theta - My_imag * cos_theta) * cos_theta;
		Qz_imag = (My_imag * cos_theta - Mz_imag * sin_theta) * sin_theta;

		// nuclear SANS cross section projected in (qz, qy)-plane
		SANSData.S_Nuc_2D_unpolarized[i] = v * (Nuc_real * Nuc_real + Nuc_imag * Nuc_imag);

		// unpolarized magnetic SANS cross section projected in (qz, qy)-plane
		SANSData.S_Mag_2D_unpolarized[i] = v * (Qx_real * Qx_real + Qx_imag * Qx_imag) \
										 + v * (Qy_real * Qy_real + Qy_imag * Qy_imag) \
										 + v * (Qz_real * Qz_real + Qz_imag * Qz_imag);

		// nuclear magnetic interference SANS cross section projected in (qz, qy)-plane
		SANSData.S_NucMag_2D[i] = 2.0 * v * Px * (Nuc_real * Qx_real + Nuc_imag * Qx_imag) \
								+ 2.0 * v * Py * (Nuc_real * Qy_real + Nuc_imag * Qy_imag) \
								+ 2.0 * v * Pz * (Nuc_real * Qz_real + Nuc_imag * Qz_imag);

		// polarized magnetic SANS cross section projected in the (qz, qy)-plane
		SANSData.S_Mag_2D_polarized[i] = v * pow(Px, 2) * (Qx_real * Qx_real + Qx_imag * Qx_imag) \
									   + v * pow(Py, 2) * (Qy_real * Qy_real + Qy_imag * Qy_imag) \
									   + v * pow(Pz, 2) * (Qz_real * Qz_real + Qz_imag * Qz_imag) \
									   + v * 2.0 * Px * Py * (Qx_real * Qy_real + Qx_imag * Qy_imag) \
									   + v * 2.0 * Px * Pz * (Qx_real * Qz_real + Qx_imag * Qz_imag) \
									   + v * 2.0 * Py * Pz * (Qy_real * Qz_real + Qy_imag * Qz_imag);

		// chiral magnetic SANS cross section in (qz, qy)-plane
		SANSData.S_Mag_2D_chiral[i] = v * 2.0 * Px * (Qy_imag * Qz_real - Qz_imag * Qy_real) \
									+ v * 2.0 * Py * (Qz_imag * Qx_real - Qx_imag * Qz_real) \
									+ v * 2.0 * Pz * (Qx_imag * Qy_real - Qy_imag * Qx_real);


		SANSData.Gxx_real[i] = v*(Mx_real * Mx_real + Mx_imag * Mx_imag);
		SANSData.Gxx_imag[i] = 0.0;

		SANSData.Gyy_real[i] = v*(My_real * My_real + My_imag * My_imag);
		SANSData.Gyy_imag[i] = 0.0;

		SANSData.Gzz_real[i] = v*(Mz_real * Mz_real + Mz_imag * Mz_imag);
		SANSData.Gzz_imag[i] = 0.0;

		SANSData.Gxy_real[i] = v*(Mx_real * My_real + Mx_imag * My_imag);
		SANSData.Gxy_imag[i] = v*(Mx_imag * My_real - Mx_real * My_imag);

		SANSData.Gyx_real[i] =  SANSData.Gxy_real[i];
		SANSData.Gyx_imag[i] = -SANSData.Gxy_imag[i];

		SANSData.Gxz_real[i] = v*(Mx_real * Mz_real + Mx_imag * Mz_imag);
		SANSData.Gxz_imag[i] = v*(Mx_imag * Mz_real - Mx_real * Mz_imag);

		SANSData.Gzx_real[i] =  SANSData.Gxz_real[i];
		SANSData.Gzx_imag[i] = -SANSData.Gxy_imag[i];

		SANSData.Gyz_real[i] = v*(My_real * Mz_real + My_imag * Mz_imag);
		SANSData.Gyz_imag[i] = v*(My_imag * Mz_real - My_real * Mz_imag);

		SANSData.Gzx_real[i] =  SANSData.Gyz_real[i];
		SANSData.Gzx_imag[i] = -SANSData.Gyz_imag[i];

	}
}
