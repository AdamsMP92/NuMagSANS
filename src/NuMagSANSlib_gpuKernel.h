// File         : NuMagSANSlib_gpuKernel.h
// Author       : Michael Philipp ADAMS, M.Sc.
// Company      : University of Luxembourg
// Department   : Department of Physics and Materials Sciences
// Group        : NanoMagnetism Group
// Group Leader : Prof. Andreas Michels
// Version      : 22 November 2024
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
             SANSData.S_Mag_1D_spin_flip[i] += SANSData.S_Mag_2D_spin_flip[j + i*N_theta] \
                                             + SANSData.S_Mag_2D_spin_flip[j + i*N_theta + 1];
             SANSData.S_Mag_1D_chiral[i] += abs(SANSData.S_Mag_2D_chiral[j + i*N_theta]) \
                                          + abs(SANSData.S_Mag_2D_chiral[j + i*N_theta + 1]);
          }
         SANSData.S_Nuc_1D_unpolarized[i] = SANSData.S_Nuc_1D_unpolarized[i]/(4.0*M_PI)*dtheta;
         SANSData.S_Mag_1D_unpolarized[i] = SANSData.S_Mag_1D_unpolarized[i]/(4.0*M_PI)*dtheta;
         SANSData.S_Mag_1D_spin_flip[i] = SANSData.S_Mag_1D_spin_flip[i]/(4.0*M_PI)*dtheta;
         SANSData.S_Mag_1D_chiral[i] = SANSData.S_Mag_1D_chiral[i]/(4.0*M_PI)*dtheta;
     }
}



// computes in the first step the correlation function c(r) and the by multiplication with r^2 the pair-distance function
 // here we take into account the limit of sin(x)/x at x-> 0 and so the singularity is fixed
 __global__
 void DistributionFunctions(ScatteringData SANSData){
  
       int i = blockIdx.x * blockDim.x + threadIdx.x;

       unsigned int N_r = *SANSData.N_r;
       unsigned int N_q = *SANSData.N_q;
       float dq = *SANSData.dq;
       float qr1 = 0.0;
       float qr2 = 0.0;
       float po1 = 0.0;
       float po2 = 0.0;
       float s1 = 0.0;
       float s2 = 0.0;
  
       if(i < N_r){
           if(i != 0){
               qr1 = SANSData.q_1D[1] * SANSData.r_1D[i];
               po1 = pow(SANSData.q_1D[1], 2);
               s1 = sin(qr1)/qr1 * po1;
               SANSData.c_Nuc_unpolarized[i] += SANSData.S_Nuc_1D_unpolarized[1] * s1;
               for(int j=1; j<N_q-1; j++){
  
                   qr1 = SANSData.q_1D[j] * SANSData.r_1D[i];
                   po1 = pow(SANSData.q_1D[j], 2);
                   s1 = sin(qr1)/qr1 * po1;
	               qr2 = SANSData.q_1D[j+1] * SANSData.r_1D[i];
	 	           po2 = pow(SANSData.q_1D[j+1], 2);
        	       s2 = sin(qr2)/qr2 * po2;

             	   SANSData.c_Nuc_unpolarized[i] += SANSData.S_Nuc_1D_unpolarized[j]  * s1  + SANSData.S_Nuc_1D_unpolarized[j+1] * s2;
               	   SANSData.c_Mag_unpolarized[i] += SANSData.S_Mag_1D_unpolarized[j]  * s1  + SANSData.S_Mag_1D_unpolarized[j+1] * s2;
               	   SANSData.c_Mag_spin_flip[i] += SANSData.S_Mag_1D_spin_flip[j]  * s1  + SANSData.S_Mag_1D_spin_flip[j+1] * s2;
           		}

		        SANSData.c_Nuc_unpolarized[i] = SANSData.c_Nuc_unpolarized[i]/2.0 * dq;
                SANSData.p_Nuc_unpolarized[i] = SANSData.c_Nuc_unpolarized[i] * pow(SANSData.r_1D[i], 2);

                SANSData.c_Mag_unpolarized[i] = SANSData.c_Mag_unpolarized[i]/2.0 * dq;
                SANSData.p_Mag_unpolarized[i] = SANSData.c_Mag_unpolarized[i] * pow(SANSData.r_1D[i], 2);
   
                SANSData.c_Mag_spin_flip[i] = SANSData.c_Mag_spin_flip[i]/2.0 * dq;
                SANSData.p_Mag_spin_flip[i] = SANSData.c_Mag_spin_flip[i] * pow(SANSData.r_1D[i], 2);
           }
           else{
                SANSData.c_Nuc_unpolarized[i] += SANSData.S_Nuc_1D_unpolarized[1] * pow(SANSData.q_1D[1], 2);
                SANSData.c_Mag_unpolarized[i] += SANSData.S_Mag_1D_unpolarized[1] * pow(SANSData.q_1D[1], 2);
                SANSData.c_Mag_spin_flip[i] += SANSData.S_Mag_1D_spin_flip[1] * pow(SANSData.q_1D[1], 2);
                for(int j=1; j<N_q-1; j++){
   
                    SANSData.c_Nuc_unpolarized[i] += SANSData.S_Nuc_1D_unpolarized[j]   * pow(SANSData.q_1D[j],   2) \
                                                   + SANSData.S_Nuc_1D_unpolarized[j+1] * pow(SANSData.q_1D[j+1], 2);
   
                    SANSData.c_Mag_unpolarized[i] += SANSData.S_Mag_1D_unpolarized[j]   * pow(SANSData.q_1D[j],   2) \
                                                   + SANSData.S_Mag_1D_unpolarized[j+1] * pow(SANSData.q_1D[j+1], 2);

                    SANSData.c_Mag_spin_flip[i] += SANSData.S_Mag_1D_spin_flip[j]   * pow(SANSData.q_1D[j],   2) \
                                                 + SANSData.S_Mag_1D_spin_flip[j+1] * pow(SANSData.q_1D[j+1], 2);
   
               }
   
               SANSData.c_Nuc_unpolarized[i] = SANSData.c_Nuc_unpolarized[i]/2.0 * dq;
               SANSData.p_Nuc_unpolarized[i] = SANSData.c_Nuc_unpolarized[i] * pow(SANSData.r_1D[i], 2);
  
               SANSData.c_Mag_unpolarized[i] = SANSData.c_Mag_unpolarized[i]/2.0 * dq;
               SANSData.p_Mag_unpolarized[i] = SANSData.c_Mag_unpolarized[i] * pow(SANSData.r_1D[i], 2);
  
               SANSData.c_Mag_spin_flip[i] = SANSData.c_Mag_spin_flip[i]/2.0 * dq;
               SANSData.p_Mag_spin_flip[i] = SANSData.c_Mag_spin_flip[i] * pow(SANSData.r_1D[i], 2);
  
           }
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
            SANSData.Corr_Mag_2D_spin_flip[i] += v * SANSData.S_Mag_2D_spin_flip[k] * c;
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
	unsigned long int N = *MagData.N;

	float v = 1.0/((float)  (*MagData.K)) * pow(1.0/((float) (*MagData.N)), 2); // pre factor
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
 
	if(i < L){
		for(int k=0; k< (*MagData.K); k++){

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

        	for(int l=0; l<N; l++){
             //X = RotMat[0] * x[l+k*N] + RotMat[3] * y[l+k*N] + RotMat[6] * z[l+k*N];
            	Y = MagData.RotMat[1] * MagData.x[l+k*N]  + MagData.RotMat[4] * MagData.y[l+k*N] + MagData.RotMat[7] * MagData.z[l+k*N];
            	Z = MagData.RotMat[2] * MagData.x[l+k*N]  + MagData.RotMat[5] * MagData.y[l+k*N] + MagData.RotMat[8] * MagData.z[l+k*N];

            	mx_real += MagData.mx[l+k*N] * cos(Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i]);
            	mx_imag -= MagData.mx[l+k*N] * sin(Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i]);
            	my_real += MagData.my[l+k*N] * cos(Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i]);
            	my_imag -= MagData.my[l+k*N] * sin(Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i]);
            	mz_real += MagData.mz[l+k*N] * cos(Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i]);
            	mz_imag -= MagData.mz[l+k*N] * sin(Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i]);
			}

			Mx_real = MagData.RotMat[0] * mx_real + MagData.RotMat[3] * my_real + MagData.RotMat[6] * mz_real;
			My_real = MagData.RotMat[1] * mx_real + MagData.RotMat[4] * my_real + MagData.RotMat[7] * mz_real;
			Mz_real = MagData.RotMat[2] * mx_real + MagData.RotMat[5] * my_real + MagData.RotMat[8] * mz_real;

			Mx_imag = MagData.RotMat[0] * mx_imag + MagData.RotMat[3] * my_imag + MagData.RotMat[6] * mz_imag;
			My_imag = MagData.RotMat[1] * mx_imag + MagData.RotMat[4] * my_imag + MagData.RotMat[7] * mz_imag;
			Mz_imag = MagData.RotMat[2] * mx_imag + MagData.RotMat[5] * my_imag + MagData.RotMat[8] * mz_imag;

            Qx_real = (-Mx_real);
            Qy_real = (Mz_real * sin(SANSData.theta_2D[i]) - My_real * cos(SANSData.theta_2D[i])) * cos(SANSData.theta_2D[i]);
            Qz_real = (My_real * cos(SANSData.theta_2D[i]) - Mz_real * sin(SANSData.theta_2D[i])) * sin(SANSData.theta_2D[i]);

            Qx_imag = (-Mx_imag);
            Qy_imag = (Mz_imag * sin(SANSData.theta_2D[i]) - My_imag * cos(SANSData.theta_2D[i])) * cos(SANSData.theta_2D[i]);
            Qz_imag = (My_imag * cos(SANSData.theta_2D[i]) - Mz_imag * sin(SANSData.theta_2D[i])) * sin(SANSData.theta_2D[i]);


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

	// spin-flip magnetic SANS cross section projected in (qz, qy)-plane
    SANSData.S_Mag_2D_spin_flip[i] = SANSData.S_Mag_2D_unpolarized[i] - SANSData.S_Mag_2D_polarized[i];

	// pm-spin-flip magnetic SANS cross section projected in (qz, qy)-plane
	SANSData.S_Mag_2D_spin_flip_pm[i] = SANSData.S_Mag_2D_spin_flip[i] + SANSData.S_Mag_2D_chiral[i];

    //mp-spin-flip magnetic SANS cross section projected in (qz, qy)-plane
	SANSData.S_Mag_2D_spin_flip_mp[i] = SANSData.S_Mag_2D_spin_flip[i] - SANSData.S_Mag_2D_chiral[i];

    // non-spin-flip magnetic SANS cross section projected in (qz, qy)-plane
    SANSData.S_Mag_2D_non_spin_flip_pp[i] = SANSData.S_Nuc_2D_unpolarized[i] + SANSData.S_NucMag_2D[i] + SANSData.S_Mag_2D_polarized[i];
    SANSData.S_Mag_2D_non_spin_flip_mm[i] = SANSData.S_Nuc_2D_unpolarized[i] - SANSData.S_NucMag_2D[i] + SANSData.S_Mag_2D_polarized[i];


    // sanspol cross sections projected in (qz, qy)-plane
    SANSData.S_Mag_2D_sanspol_p[i] = SANSData.S_Mag_2D_non_spin_flip_pp[i] + SANSData.S_Mag_2D_spin_flip_pm[i];
    SANSData.S_Mag_2D_sanspol_m[i] = SANSData.S_Mag_2D_non_spin_flip_mm[i] + SANSData.S_Mag_2D_spin_flip_mp[i];

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
	unsigned long int N = *NucData.N;

	float v = 1.0/((float)  (*NucData.K)) * pow(1.0/((float) (*NucData.N)), 2); // pre factor

	int i = blockIdx.x * blockDim.x + threadIdx.x;

	float Nuc_real = 0.0;
	float Nuc_imag = 0.0;
	float Y = 0.0;
	float Z = 0.0;
 
	if(i < L){
		for(int k=0; k< (*NucData.K); k++){

        	Nuc_real = 0.0;
        	Nuc_imag = 0.0;

        	for(int l=0; l<N; l++){
             //X = RotMat[0] * x[l+k*N] + RotMat[3] * y[l+k*N] + RotMat[6] * z[l+k*N];
            	Y = NucData.RotMat[1] * NucData.x[l+k*N]  + NucData.RotMat[4] * NucData.y[l+k*N] + NucData.RotMat[7] * NucData.z[l+k*N];
            	Z = NucData.RotMat[2] * NucData.x[l+k*N]  + NucData.RotMat[5] * NucData.y[l+k*N] + NucData.RotMat[8] * NucData.z[l+k*N];

				Nuc_real += NucData.Nuc[l+k*N] * cos(Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i]);
				Nuc_imag -= NucData.Nuc[l+k*N] * sin(Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i]);
            	
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
	unsigned long int N = *MagData.N;
	unsigned long int K = *MagData.K;

	float v = (1.0/((float) K)) * pow(1.0/((float) N), 2); // pre factor
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

	float Nuc_real = 0.0;
	float Nuc_imag = 0.0;
   // float X = 0.0;
	float Y = 0.0;
	float Z = 0.0;

	float Psi = 0.0;

	float cos_val = 0.0;
	float sin_val = 0.0;
 
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

        	Nuc_real = 0.0;
        	Nuc_imag = 0.0;

        	for(int l=0; l<N; l++){
				// atomic position composition
				//X = MagData.RotMat[0] * MagData.x[l+k*N] \
                //  + MagData.RotMat[3] * MagData.y[l+k*N] \
                //  + MagData.RotMat[6] * MagData.z[l+k*N];
            	Y = MagData.RotMat[1] * MagData.x[l+k*N] \
            	  + MagData.RotMat[4] * MagData.y[l+k*N] \
            	  + MagData.RotMat[7] * MagData.z[l+k*N];
            	Z = MagData.RotMat[2] * MagData.x[l+k*N] \
            	  + MagData.RotMat[5] * MagData.y[l+k*N] \
            	  + MagData.RotMat[8] * MagData.z[l+k*N];

				// phase function
				Psi = Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i];

				// cosine and sine values
				cos_val = cos(Psi);
				sin_val = sin(Psi);

				// cosine and sine summations
				Nuc_real += NucData.Nuc[l+k*N] * cos_val;
				Nuc_imag -= NucData.Nuc[l+k*N] * sin_val;

            	mx_real += MagData.mx[l+k*N] * cos_val;
            	mx_imag -= MagData.mx[l+k*N] * sin_val;
            	my_real += MagData.my[l+k*N] * cos_val;
            	my_imag -= MagData.my[l+k*N] * sin_val;
            	mz_real += MagData.mz[l+k*N] * cos_val;
            	mz_imag -= MagData.mz[l+k*N] * sin_val;




				/*
                X = RotMat[0] * x[l+k*N] + RotMat[3] * y[l+k*N] + RotMat[6] * z[l+k*N];
				Y = MagData.RotMat[1] * MagData.x[l+k*N]  + MagData.RotMat[4] * MagData.y[l+k*N] + MagData.RotMat[7] * MagData.z[l+k*N];
            	Z = MagData.RotMat[2] * MagData.x[l+k*N]  + MagData.RotMat[5] * MagData.y[l+k*N] + MagData.RotMat[8] * MagData.z[l+k*N];

				Nuc_real += NucData.Nuc[l+k*N] * cos(Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i]);
				Nuc_imag -= NucData.Nuc[l+k*N] * sin(Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i]);

            	mx_real += MagData.mx[l+k*N] * cos(Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i]);
            	mx_imag -= MagData.mx[l+k*N] * sin(Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i]);
            	my_real += MagData.my[l+k*N] * cos(Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i]);
            	my_imag -= MagData.my[l+k*N] * sin(Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i]);
            	mz_real += MagData.mz[l+k*N] * cos(Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i]);
            	mz_imag -= MagData.mz[l+k*N] * sin(Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i]);
            	*/

			}

			// rotations of the magnetization fourier components real part
			Mx_real = MagData.RotMat[0] * mx_real + MagData.RotMat[3] * my_real + MagData.RotMat[6] * mz_real;
			My_real = MagData.RotMat[1] * mx_real + MagData.RotMat[4] * my_real + MagData.RotMat[7] * mz_real;
			Mz_real = MagData.RotMat[2] * mx_real + MagData.RotMat[5] * my_real + MagData.RotMat[8] * mz_real;

			// rotations of the magnetization fourier components imaginary part
			Mx_imag = MagData.RotMat[0] * mx_imag + MagData.RotMat[3] * my_imag + MagData.RotMat[6] * mz_imag;
			My_imag = MagData.RotMat[1] * mx_imag + MagData.RotMat[4] * my_imag + MagData.RotMat[7] * mz_imag;
			Mz_imag = MagData.RotMat[2] * mx_imag + MagData.RotMat[5] * my_imag + MagData.RotMat[8] * mz_imag;

			// real-parts of the Halpern-Johnson vector
			Qx_real = (-Mx_real);
            Qy_real = (Mz_real * sin(SANSData.theta_2D[i]) - My_real * cos(SANSData.theta_2D[i])) * cos(SANSData.theta_2D[i]);
            Qz_real = (My_real * cos(SANSData.theta_2D[i]) - Mz_real * sin(SANSData.theta_2D[i])) * sin(SANSData.theta_2D[i]);

			// imaginary-parts of the Halpern-Johnson vector
            Qx_imag = (-Mx_imag);
            Qy_imag = (Mz_imag * sin(SANSData.theta_2D[i]) - My_imag * cos(SANSData.theta_2D[i])) * cos(SANSData.theta_2D[i]);
            Qz_imag = (My_imag * cos(SANSData.theta_2D[i]) - Mz_imag * sin(SANSData.theta_2D[i])) * sin(SANSData.theta_2D[i]);



			// nuclear SANS cross section projected in (qz, qy)-plane
			SANSData.S_Nuc_2D_unpolarized[i] += v * (Nuc_real * Nuc_real + Nuc_imag * Nuc_imag);

			// unpolarized magnetic SANS cross section projected in (qz, qy)-plane
			SANSData.S_Mag_2D_unpolarized[i] += v * (Qx_real * Qx_real + Qx_imag * Qx_imag) \
											  + v * (Qy_real * Qy_real + Qy_imag * Qy_imag) \
											  + v * (Qz_real * Qz_real + Qz_imag * Qz_imag);

			// nuclear magnetic interference SANS cross section projected in (qz, qy)-plane
			SANSData.S_NucMag_2D[i] += 2.0 * v * Px * (Nuc_real * Qx_real + Nuc_imag * Qx_imag) \
									 + 2.0 * v * Py * (Nuc_real * Qy_real + Nuc_imag * Qy_imag) \
									 + 2.0 * v * Pz * (Nuc_real * Qz_real + Nuc_imag * Qz_imag);

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

    // spin-flip magnetic SANS cross section projected in (qz, qy)-plane
    SANSData.S_Mag_2D_spin_flip[i] = SANSData.S_Mag_2D_unpolarized[i] - SANSData.S_Mag_2D_polarized[i];

	// pm-spin-flip magnetic SANS cross section projected in (qz, qy)-plane
	SANSData.S_Mag_2D_spin_flip_pm[i] = SANSData.S_Mag_2D_spin_flip[i] + SANSData.S_Mag_2D_chiral[i];

    //mp-spin-flip magnetic SANS cross section projected in (qz, qy)-plane
	SANSData.S_Mag_2D_spin_flip_mp[i] = SANSData.S_Mag_2D_spin_flip[i] - SANSData.S_Mag_2D_chiral[i];

    // non-spin-flip magnetic SANS cross section projected in (qz, qy)-plane
    SANSData.S_Mag_2D_non_spin_flip_pp[i] = SANSData.S_Nuc_2D_unpolarized[i] + SANSData.S_NucMag_2D[i] + SANSData.S_Mag_2D_polarized[i];
    SANSData.S_Mag_2D_non_spin_flip_mm[i] = SANSData.S_Nuc_2D_unpolarized[i] - SANSData.S_NucMag_2D[i] + SANSData.S_Mag_2D_polarized[i];


    // sanspol cross sections projected in (qz, qy)-plane
    SANSData.S_Mag_2D_sanspol_p[i] = SANSData.S_Mag_2D_non_spin_flip_pp[i] + SANSData.S_Mag_2D_spin_flip_pm[i];
    SANSData.S_Mag_2D_sanspol_m[i] = SANSData.S_Mag_2D_non_spin_flip_mm[i] + SANSData.S_Mag_2D_spin_flip_mp[i];


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
	unsigned long int N = *MagData.N;
	unsigned long int K = *MagData.K;

	float v = (1.0/((float) K)) * pow(1.0/((float) N), 2); // pre factor
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

	float Nuc_real = 0.0;
	float Nuc_imag = 0.0;
   // float X = 0.0;
	float Y = 0.0;
	float Z = 0.0;

	float Psi = 0.0;

	float cos_val = 0.0;
	float sin_val = 0.0;


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

        	Nuc_real = 0.0;
        	Nuc_imag = 0.0;

        	for(int l=0; l<N; l++){

				// atomic position composition
				//X = MagData.RotMat[0] * (MagData.x[l+k*N] + StructData.x[k]) \
                //  + MagData.RotMat[3] * (MagData.y[l+k*N] + StructData.y[k]) \
				/ / + MagData.RotMat[6] * (MagData.z[l+k*N] + StructData.z[k]);
            	Y = MagData.RotMat[1] * (MagData.x[l+k*N] + StructData.x[k]) \
            	  + MagData.RotMat[4] * (MagData.y[l+k*N] + StructData.y[k]) \
            	  + MagData.RotMat[7] * (MagData.z[l+k*N] + StructData.z[k]);
            	Z = MagData.RotMat[2] * (MagData.x[l+k*N] + StructData.x[k]) \
            	  + MagData.RotMat[5] * (MagData.y[l+k*N] + StructData.y[k]) \
            	  + MagData.RotMat[8] * (MagData.z[l+k*N] + StructData.z[k]);

				// phase function
				Psi = Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i];

				// cosine and sine values
				cos_val = cos(Psi);
				sin_val = sin(Psi);

				// cosine and sine summations
				Nuc_real += NucData.Nuc[l+k*N] * cos_val;
				Nuc_imag -= NucData.Nuc[l+k*N] * sin_val;

            	mx_real += MagData.mx[l+k*N] * cos_val;
            	mx_imag -= MagData.mx[l+k*N] * sin_val;
            	my_real += MagData.my[l+k*N] * cos_val;
            	my_imag -= MagData.my[l+k*N] * sin_val;
            	mz_real += MagData.mz[l+k*N] * cos_val;
            	mz_imag -= MagData.mz[l+k*N] * sin_val;
			}

			// rotations of the magnetization fourier components real part
			Mx_real += MagData.RotMat[0] * mx_real + MagData.RotMat[3] * my_real + MagData.RotMat[6] * mz_real;
			My_real += MagData.RotMat[1] * mx_real + MagData.RotMat[4] * my_real + MagData.RotMat[7] * mz_real;
			Mz_real += MagData.RotMat[2] * mx_real + MagData.RotMat[5] * my_real + MagData.RotMat[8] * mz_real;

			// rotations of the magnetization fourier components imaginary part
			Mx_imag += MagData.RotMat[0] * mx_imag + MagData.RotMat[3] * my_imag + MagData.RotMat[6] * mz_imag;
			My_imag += MagData.RotMat[1] * mx_imag + MagData.RotMat[4] * my_imag + MagData.RotMat[7] * mz_imag;
			Mz_imag += MagData.RotMat[2] * mx_imag + MagData.RotMat[5] * my_imag + MagData.RotMat[8] * mz_imag;

		}

		// real-parts of the Halpern-Johnson vector
		Qx_real = (-Mx_real);
		Qy_real = (Mz_real * sin(SANSData.theta_2D[i]) - My_real * cos(SANSData.theta_2D[i])) * cos(SANSData.theta_2D[i]);
		Qz_real = (My_real * cos(SANSData.theta_2D[i]) - Mz_real * sin(SANSData.theta_2D[i])) * sin(SANSData.theta_2D[i]);

		// imaginary-parts of the Halpern-Johnson vector
		Qx_imag = (-Mx_imag);
		Qy_imag = (Mz_imag * sin(SANSData.theta_2D[i]) - My_imag * cos(SANSData.theta_2D[i])) * cos(SANSData.theta_2D[i]);
		Qz_imag = (My_imag * cos(SANSData.theta_2D[i]) - Mz_imag * sin(SANSData.theta_2D[i])) * sin(SANSData.theta_2D[i]);



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



		// spin-flip magnetic SANS cross section projected in (qz, qy)-plane
		SANSData.S_Mag_2D_spin_flip[i] = SANSData.S_Mag_2D_unpolarized[i] - SANSData.S_Mag_2D_polarized[i];

		// pm-spin-flip magnetic SANS cross section projected in (qz, qy)-plane
		SANSData.S_Mag_2D_spin_flip_pm[i] = SANSData.S_Mag_2D_spin_flip[i] + SANSData.S_Mag_2D_chiral[i];

		//mp-spin-flip magnetic SANS cross section projected in (qz, qy)-plane
		SANSData.S_Mag_2D_spin_flip_mp[i] = SANSData.S_Mag_2D_spin_flip[i] - SANSData.S_Mag_2D_chiral[i];

		// non-spin-flip magnetic SANS cross section projected in (qz, qy)-plane
		SANSData.S_Mag_2D_non_spin_flip_pp[i] = SANSData.S_Nuc_2D_unpolarized[i] + SANSData.S_NucMag_2D[i] + SANSData.S_Mag_2D_polarized[i];
		SANSData.S_Mag_2D_non_spin_flip_mm[i] = SANSData.S_Nuc_2D_unpolarized[i] - SANSData.S_NucMag_2D[i] + SANSData.S_Mag_2D_polarized[i];


		// sanspol cross sections projected in (qz, qy)-plane
		SANSData.S_Mag_2D_sanspol_p[i] = SANSData.S_Mag_2D_non_spin_flip_pp[i] + SANSData.S_Mag_2D_spin_flip_pm[i];
		SANSData.S_Mag_2D_sanspol_m[i] = SANSData.S_Mag_2D_non_spin_flip_mm[i] + SANSData.S_Mag_2D_spin_flip_mp[i];


	}
}











// gpu implementation for the multiplication of the Fourier transform result with the micromagnetic cuboid form factor function
__global__
void Micromagnetic_FormFactor_Multiplier(int L, \
                    float*qy, float*qz, \
                   float cuboid_cell_size_x, \
                  float cuboid_cell_size_y, \
                  float cuboid_cell_size_z, \
                   float*Gxx_real, float*Gyy_real, float*Gzz_real, float*Gxy_real, float*Gyx_real, float*Gxz_real, float*Gzx_real, float*Gyz_real, float*Gzy_real, \
                   float*Gxx_imag, float*Gyy_imag, float*Gzz_imag, float*Gxy_imag, float*Gyx_imag, float*Gxz_imag, float*Gzx_imag, float*Gyz_imag, float*Gzy_imag, \
                        float*dSigma_dOmega_M_unpolarized, \
                      float*dSigma_dOmega_M_spin_flip, \
                       float*dSigma_dOmega_M_chiral, \
                      float*dSigma_dOmega_M_spin_flip_pm, \
                       float*dSigma_dOmega_M_spin_flip_mp, \
                       float*RotMat){

    // Input information:
    // L     : number of points in Fourier space L = N_q*N_theta
    // qy    : qy-Fourier-space coordinate in units of inverse nano-meters
     // qz    : qz-Fourier-space coordinate in units of inverse nano-meters

    // output information:
    // Gxx_real: real-part of the xx-component of the Fourier-space correlation function of the magnetization
    // Gxx_imag: imaginary-part of the xx-component of the Fourier-space correlation function of the magnetization ...

    // dSigma_dOmega_M_unpolarized  : magnetic unpolarized SANS cross section
   // dSigma_dOmega_M_spin_flip    : magnetic spin-flip SANS cross section (sum of pm and mp )/2
   // dSigma_dOmega_M_chiral       : magnetic chiral SANS cross section (difference of pm and mp)/2
    // dSigma_dOmega_M_spin_flip_pm : magnetic pm spin-flip SANS cross section
    // dSigma_dOmega_M_spin_flip_mp : magnetic mp spin-flip SANS cross section


	int i = blockIdx.x * blockDim.x + threadIdx.x;
	float qx = 0.0;
	float Qx = 0.0;
	float Qy = 0.0;
	float Qz = 0.0;
	float jx = 0.0;
	float jy = 0.0;
	float jz = 0.0;
	float jxyz = 0.0;

 	if(i < L){
 
		Qx = RotMat[0] * qx + RotMat[1] * qy[i] + RotMat[2] * qz[i];
		Qy = RotMat[3] * qx + RotMat[4] * qy[i] + RotMat[5] * qz[i];
		Qz = RotMat[6] * qx + RotMat[7] * qy[i] + RotMat[8] * qz[i];

		if(Qx != 0.0){
			jx = sin(Qx*cuboid_cell_size_x/2)/(Qx*cuboid_cell_size_x/2);
		}
		else{
			jx = 1.0;
		}

		if(Qy != 0.0){
			jy = sin(Qy*cuboid_cell_size_y/2)/(Qy*cuboid_cell_size_y/2);
		}
		else{
			jy = 1.0;
		}

		if(qz[i] != 0.0){
			jz = sin(Qz*cuboid_cell_size_z/2)/(Qz*cuboid_cell_size_z/2);
		}
        else{
			jz = 1.0;
		}

		jxyz = pow(jz * jy * jx, 2);

		Gxx_real[i] = jxyz*Gxx_real[i];
		Gxx_imag[i] = jxyz*Gxx_imag[i];

		Gyy_real[i] = jxyz*Gyy_real[i];
		Gyy_imag[i] = jxyz*Gyy_imag[i];

		Gzz_real[i] = jxyz*Gzz_real[i];
		Gzz_imag[i] = jxyz*Gzz_imag[i];

		Gxy_real[i] = jxyz*Gxy_real[i];
		Gxy_imag[i] = jxyz*Gxy_imag[i];

		Gyx_real[i] = jxyz*Gyx_real[i];
		Gyx_imag[i] = jxyz*Gyx_imag[i];

		Gxz_real[i] = jxyz*Gxz_real[i];
		Gxz_imag[i] = jxyz*Gxz_imag[i];

		Gzx_real[i] = jxyz*Gzx_real[i];
		Gzx_imag[i] = jxyz*Gzx_imag[i];

		Gyz_real[i] = jxyz*Gyz_real[i];
		Gyz_imag[i] = jxyz*Gyz_imag[i];

		Gzx_real[i] = jxyz*Gzx_real[i];
		Gzx_imag[i] = jxyz*Gzx_imag[i];

		// unpolarized magnetic SANS cross section projected in (qz, qy)-plane
		dSigma_dOmega_M_unpolarized[i] = jxyz*dSigma_dOmega_M_unpolarized[i];

		// chiral magnetic SANS cross section in (qz, qy)-plane
		dSigma_dOmega_M_chiral[i] = jxyz*dSigma_dOmega_M_chiral[i];

		// spin-flip magnetic SANS cross section projected in (qz, qy)-plane
		dSigma_dOmega_M_spin_flip[i] = jxyz*dSigma_dOmega_M_spin_flip[i];

		// pm-spin-flip magnetic SANS cross section projected in (qz, qy)-plane
		dSigma_dOmega_M_spin_flip_pm[i] = jxyz*dSigma_dOmega_M_spin_flip_pm[i];

		//mp-spin-flip magnetic SANS cross section projected in (qz, qy)-plane
		dSigma_dOmega_M_spin_flip_mp[i] = jxyz*dSigma_dOmega_M_spin_flip_mp[i];

	}
}

