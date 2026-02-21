// File         : NuMagSANSlib_SANSData.h
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

using namespace std;

struct ScatteringData {

	unsigned int *N_q;		// number of q-values
	unsigned int *N_theta;	// number of theta-values
	unsigned int *N_r;		// number of r-values
	unsigned int *N_alpha;	// number of alpha-values

	float *Polarization; // polarization vector [Px, Py, Pz]

	float *q_max;	// maximum q-value
	float *r_max;	// maximum r-value
	
	float *dq;		// step-size of q
	float *dtheta;	// step-size of theta
	float *dr;		// step-size of r
	float *dalpha;	// step-size of alpha

	float *qy_2D;  // 2D qy-scattering vector
	float *qz_2D;  // 2D qz-scattering vector
	float *q_2D;   // 2D q-scattering vector
	float *theta_2D; // 2D theta-scattering vector

	float *ry_2D;		// 2D ry-correlation space coordinate
	float *rz_2D;   	// 2D rz-correlation space coordinate
	float *r_2D;		// 2D r-correlation space coordinate
	float *alpha_2D;	// 2D angle-correlation space coordinate

	float *q_1D;	// 1D scattering vector
	float *r_1D;	// 1D correlation space coordinate

	float *Gxx_real; // real-part xx-magnetic Fourier correlation function
	float *Gyy_real; // real-part yy-magnetic Fourier correlation function
	float *Gzz_real; // real part zz-magnetic Fourier correlation function
	float *Gxy_real; // real part xy-magnetic Fourier correlation function
	float *Gyx_real; // real part yx-magnetic Fourier correlation function
	float *Gxz_real; // real part xz-magnetic Fourier correlation function
	float *Gzx_real; // real part zx-magnetic Fourier correlation function
	float *Gyz_real; // real part yz-magnetic Fourier correlation function
	float *Gzy_real; // real part zy-magnetic Fourier correlation function

	float *Gxx_imag; // real-part xx-magnetic Fourier correlation function
	float *Gyy_imag; // real-part yy-magnetic Fourier correlation function
	float *Gzz_imag; // real-part zz-magnetic Fourier correlation function
	float *Gxy_imag; // real-part xy-magnetic Fourier correlation function
	float *Gyx_imag; // real-part yx-magnetic Fourier correlation function
	float *Gxz_imag; // real-part xz-magnetic Fourier correlation function
	float *Gzx_imag; // real-part zx-magnetic Fourier correlation function
	float *Gyz_imag; // real-part yz-magnetic Fourier correlation function
	float *Gzy_imag; // real-part zy-magnetic Fourier correlation function

	float *S_Nuc_2D_unpolarized;		// nuclear SANS cross section
	float *S_Mag_2D_unpolarized;		// unpolarized magnetic SANS cross section
	float *S_Mag_2D_polarized;			// polarized magnetic SANS cross section
	float *S_NucMag_2D;					// nuclear-magnetic interference SANS cross section
	float *S_Mag_2D_spin_flip;			// spin-flip magnetic SANS cross section
	float *S_Mag_2D_chiral;				// chiral magnetic SANS cross section
	float *S_Mag_2D_spin_flip_pm;		// pm-spin-flip magnetic SANS cross section
	float *S_Mag_2D_spin_flip_mp;		// mp-spin-flip magnetic SANS cross section
	float *S_Mag_2D_non_spin_flip_pp;	// pp-non-spin-flip magnetic SANS cross section
	float *S_Mag_2D_non_spin_flip_mm;	// mm-non-spin-flip magnetic SANS cross section
	float *S_Mag_2D_sanspol_p;			// p-sanspol magnetic SANS cross section
	float *S_Mag_2D_sanspol_m;			// m-sanspol magnetic SANS cross section

	float *Corr_Nuc_2D_unpolarized;
	float *Corr_Mag_2D_unpolarized;
	float *Corr_Mag_2D_polarized;
	float *Corr_NucMag_2D;
	float *Corr_Mag_2D_spin_flip;
	float *Corr_Mag_2D_chiral;
	float *Corr_Mag_2D_spin_flip_pm;
	float *Corr_Mag_2D_spin_flip_mp;
	float *Corr_Mag_2D_non_spin_flip_pp;
	float *Corr_Mag_2D_non_spin_flip_mm;
	float *Corr_Mag_2D_sanspol_p;
	float *Corr_Mag_2D_sanspol_m;

	float *S_Nuc_1D_unpolarized;
	float *S_Mag_1D_unpolarized;
	float *S_Mag_1D_polarized;
	float *S_NucMag_1D;
	float *S_Mag_1D_chiral;
	float *S_Mag_1D_spin_flip;
	float *S_Mag_1D_spin_flip_pm;
	float *S_Mag_1D_spin_flip_mp;
	float *S_Mag_1D_non_spin_flip_pp;
	float *S_Mag_1D_non_spin_flip_mm;
	float *S_Mag_1D_sanspol_p;
	float *S_Mag_1D_sanspol_m;

	float *p_Nuc_unpolarized;
	float *p_Mag_unpolarized;
	float *p_Mag_polarized;
	float *p_NucMag;
	float *p_Mag_chiral;
	float *p_Mag_spin_flip;
	float *p_Mag_spin_flip_pm;
	float *p_Mag_spin_flip_mp;
	float *p_Mag_non_spin_flip_pp;
	float *p_Mag_non_spin_flip_mm;
	float *p_Mag_sanspol_p;
	float *p_Mag_sanspol_m;

	float *c_Nuc_unpolarized;
	float *c_Mag_unpolarized;
	float *c_Mag_polarized;
	float *c_NucMag;
	float *c_Mag_chiral;
	float *c_Mag_spin_flip;
	float *c_Mag_spin_flip_pm;
	float *c_Mag_spin_flip_mp;
	float *c_Mag_non_spin_flip_pp;
	float *c_Mag_non_spin_flip_mm;
	float *c_Mag_sanspol_p;
	float *c_Mag_sanspol_m;
	
};

// #################################################################################################################################################
struct ScalingFactors{

	float Nuc_SANS_SF;		// physical scaling factor for the nuclear SANS cross section
	float Mag_SANS_SF;		// physical scaling factor for the magnetic SANS cross section
	float NucMag_SANS_SF;	// physical scaling factor for the nuclear-magnetic SANS cross section
	float CorrelationMatrix_scaling_factor;	// physical scaling factor for the Correlation matrizes in fourier space
//	float correlation_function_1D_scaling_factor;
	//float pair_distance_distribution_1D_scaling_factor;
//	float correlation_function_2D_scaling_factor;
	float CorF_Mag_2D_SF;	// physical scaling factor for the magnetic correlation functions 2D

	float Nuc_Corr_SF; 
	float Mag_Corr_SF; 
	float NucMag_Corr_SF;
	
	float Nuc_PairDist_SF;
	float Mag_PairDist_SF;
	float NucMag_PairDist_SF;
	
};


// #################################################################################################################################################
void init_ScalingFactors(ScalingFactors* ScalFactors,\
					     InputFileData* InputData,\
					     MagnetizationData* MagData, \
					     NuclearData* NucData, \
					     ScatteringData* SANSData){

		// physical constants and scaling factors
		float mu_B = 9.2740100783*1e-24;    // Bohr Magneton in units of Am^2
		float b_H = 2.91 * 1e8;             // Magnetic Scattering length in SI units
		float n_B = 1e-15; 					// nuclear scattering length scaling factor in units of femto meters

		// Total Scattering Volume in m^3
		float V = InputData->Scattering_Volume_V;

		// Micromagnetic properties
		float sx = InputData->cuboid_cell_size_x; // cell size x in nanometers
		float sy = InputData->cuboid_cell_size_y; // cell size y in nanometers
		float sz = InputData->cuboid_cell_size_z; // cell size z in nanometers
		float Ms = InputData->cell_magnetization; // saturation magnetization in units of A/m
		float mag_moment = Ms * sx * sy * sz * 1e-27;	// cell magnetic moment

		// nuclear properties
		float cell_nuclear_sld = InputData->cell_nuclear_sld;		// nuclear scattering length density in units of 1/m^2
		float n_length = cell_nuclear_sld * sx * sy * sz * 1e-27;	// nuclear scattering length in units of

		// Total Number of atoms
		unsigned long int W;
		unsigned long int N_avg;
		if(InputData->MagData_activate_flag){
			W = *MagData->TotalAtomNumber;
			N_avg = *MagData->N_avg;
		}else{
			W = *NucData->TotalAtomNumber;
			N_avg = *NucData->N_avg;
		}
	   


	if(InputData->Fourier_Approach == "atomistic"){

		ScalFactors->Nuc_SANS_SF = ((float) W) * ((float) N_avg) * pow(n_B, 2)/V * 1e-2;
		ScalFactors->Mag_SANS_SF = ((float) W) * ((float) N_avg) * pow(mu_B * b_H, 2)/V * 1e-2;  // scaling factor for the SANS cross sections
		ScalFactors->NucMag_SANS_SF = ((float) W) * ((float) N_avg) * (n_B * mu_B * b_H)/V * 1e-2;
		ScalFactors->CorrelationMatrix_scaling_factor = pow(((float) N_avg) * mu_B, 2)/(8.0*pow(M_PI, 3));   // scaling factor for the Fourier correlation functions
		
		ScalFactors->Nuc_Corr_SF = ScalFactors->Nuc_SANS_SF * 1e-7;
		ScalFactors->Mag_Corr_SF = ScalFactors->Mag_SANS_SF * 1e-7;
		ScalFactors->NucMag_Corr_SF = ScalFactors->NucMag_SANS_SF * 1e-7;
		
		ScalFactors->Nuc_PairDist_SF = ScalFactors->Nuc_SANS_SF * 1e-7;
		ScalFactors->Mag_PairDist_SF = ScalFactors->Mag_SANS_SF * 1e-7;
		ScalFactors->NucMag_PairDist_SF = ScalFactors->NucMag_SANS_SF * 1e-7;
		
		ScalFactors->CorF_Mag_2D_SF = ScalFactors->Mag_SANS_SF * (*SANSData->dr) * (*SANSData->dtheta);

	}
	else if(InputData->Fourier_Approach == "micromagnetic"){

		ScalFactors->Nuc_SANS_SF = ((float) W) * ((float) N_avg) * pow(n_length, 2)/V * 1e-2;
		ScalFactors->Mag_SANS_SF = ((float) W) * ((float) N_avg) * pow(mag_moment * b_H, 2)/V * 1e-2; // scaling factor for the SANS cross sections
		ScalFactors->NucMag_SANS_SF = ((float) W) * ((float) N_avg) * (n_length * mag_moment * b_H)/V * 1e-2;
		ScalFactors->CorrelationMatrix_scaling_factor = pow(((float) N_avg) * mag_moment, 2)/(8.0*pow(M_PI, 3));     // scaling factor for the Fourier correlation functions

		ScalFactors->Nuc_Corr_SF = ScalFactors->Nuc_SANS_SF * 1e-7;
		ScalFactors->Mag_Corr_SF = ScalFactors->Mag_SANS_SF * 1e-7;
		ScalFactors->NucMag_Corr_SF = ScalFactors->NucMag_SANS_SF * 1e-7;
		
		ScalFactors->Nuc_PairDist_SF = ScalFactors->Nuc_SANS_SF * 1e-7;
		ScalFactors->Mag_PairDist_SF = ScalFactors->Mag_SANS_SF * 1e-7;       // scaling factor for the pair-distance distribution function
		ScalFactors->NucMag_PairDist_SF = ScalFactors->NucMag_SANS_SF * 1e-7;
		
		ScalFactors->CorF_Mag_2D_SF = ScalFactors->Mag_SANS_SF * (*SANSData->dr) * (*SANSData->dtheta);

	}



	LogSystem::write("magnetic SANS scaling factor: " + std::to_string(ScalFactors->Mag_SANS_SF));
	LogSystem::write("correlation Matrix scaling factor: " + std::to_string(ScalFactors->CorrelationMatrix_scaling_factor));

	LogSystem::write("nuclear correlation function scaling factor: " + std::to_string(ScalFactors->Nuc_Corr_SF));
	LogSystem::write("magnetic correlation function scaling factor: " + std::to_string(ScalFactors->Mag_Corr_SF));
	LogSystem::write("nuclear correlation function distance scaling factor: " + std::to_string(ScalFactors->NucMag_Corr_SF));

	LogSystem::write("nuclear pair distance scaling factor: " + std::to_string(ScalFactors->Nuc_PairDist_SF));
	LogSystem::write("magnetic pair distance scaling factor: " + std::to_string(ScalFactors->Mag_PairDist_SF));
	LogSystem::write("nuclear magnetic pair distance scaling factor: " + std::to_string(ScalFactors->NucMag_PairDist_SF));

	LogSystem::write("correlation function 2D scaling factor: " + std::to_string(ScalFactors->CorF_Mag_2D_SF));
	LogSystem::write("");

}

// ###########################################################################################################################################
void scale_ScatteringData(ScalingFactors* ScalFactors, \
                          ScatteringData* SANSData, \
                          InputFileData* InputData){

	LogSystem::write("");
	LogSystem::write("scale scattering data...");
	LogSystem::write("");

//	unsigned int L = (InputData->N_q) * (InputData->N_theta);

	for(unsigned long int l=0; l < (InputData->N_q) * (InputData->N_theta); l++){
	
		SANSData->Gxx_real[l] = SANSData->Gxx_real[l] * ScalFactors->CorrelationMatrix_scaling_factor;
		SANSData->Gyy_real[l] = SANSData->Gyy_real[l] * ScalFactors->CorrelationMatrix_scaling_factor;
		SANSData->Gzz_real[l] = SANSData->Gzz_real[l] * ScalFactors->CorrelationMatrix_scaling_factor;
		SANSData->Gxy_real[l] = SANSData->Gxy_real[l] * ScalFactors->CorrelationMatrix_scaling_factor;
		SANSData->Gyx_real[l] = SANSData->Gyx_real[l] * ScalFactors->CorrelationMatrix_scaling_factor;
		SANSData->Gxz_real[l] = SANSData->Gxz_real[l] * ScalFactors->CorrelationMatrix_scaling_factor;
		SANSData->Gzx_real[l] = SANSData->Gzx_real[l] * ScalFactors->CorrelationMatrix_scaling_factor;
		SANSData->Gyz_real[l] = SANSData->Gyz_real[l] * ScalFactors->CorrelationMatrix_scaling_factor;
		SANSData->Gzy_real[l] = SANSData->Gzy_real[l] * ScalFactors->CorrelationMatrix_scaling_factor;

		SANSData->Gxx_imag[l] = SANSData->Gxx_imag[l] * ScalFactors->CorrelationMatrix_scaling_factor;
		SANSData->Gyy_imag[l] = SANSData->Gyy_imag[l] * ScalFactors->CorrelationMatrix_scaling_factor;
		SANSData->Gzz_imag[l] = SANSData->Gzz_imag[l] * ScalFactors->CorrelationMatrix_scaling_factor;
		SANSData->Gxy_imag[l] = SANSData->Gxy_imag[l] * ScalFactors->CorrelationMatrix_scaling_factor;
		SANSData->Gyx_imag[l] = SANSData->Gyx_imag[l] * ScalFactors->CorrelationMatrix_scaling_factor;
		SANSData->Gxz_imag[l] = SANSData->Gxz_imag[l] * ScalFactors->CorrelationMatrix_scaling_factor;
		SANSData->Gzx_imag[l] = SANSData->Gzx_imag[l] * ScalFactors->CorrelationMatrix_scaling_factor;
		SANSData->Gyz_imag[l] = SANSData->Gyz_imag[l] * ScalFactors->CorrelationMatrix_scaling_factor;
		SANSData->Gzy_imag[l] = SANSData->Gzy_imag[l] * ScalFactors->CorrelationMatrix_scaling_factor;

		// basic SANS cross sections
		SANSData->S_Nuc_2D_unpolarized[l] = SANSData->S_Nuc_2D_unpolarized[l] * ScalFactors->Nuc_SANS_SF;
		SANSData->S_Mag_2D_unpolarized[l] = SANSData->S_Mag_2D_unpolarized[l] * ScalFactors->Mag_SANS_SF;
		SANSData->S_Mag_2D_polarized[l] = SANSData->S_Mag_2D_polarized[l] * ScalFactors->Mag_SANS_SF;;
		SANSData->S_NucMag_2D[l] = SANSData->S_NucMag_2D[l] * ScalFactors->NucMag_SANS_SF;
		SANSData->S_Mag_2D_chiral[l] = SANSData->S_Mag_2D_chiral[l] * ScalFactors->Mag_SANS_SF;

		// composed SANS cross sections
		SANSData->S_Mag_2D_spin_flip[l] = SANSData->S_Mag_2D_unpolarized[l] - SANSData->S_Mag_2D_polarized[l];
		SANSData->S_Mag_2D_spin_flip_pm[l] = SANSData->S_Mag_2D_spin_flip[l] + SANSData->S_Mag_2D_chiral[l];
		SANSData->S_Mag_2D_spin_flip_mp[l] = SANSData->S_Mag_2D_spin_flip[l] - SANSData->S_Mag_2D_chiral[l];
		SANSData->S_Mag_2D_non_spin_flip_pp[l] = SANSData->S_Nuc_2D_unpolarized[l] + SANSData->S_NucMag_2D[l] + SANSData->S_Mag_2D_polarized[l];
		SANSData->S_Mag_2D_non_spin_flip_mm[l] = SANSData->S_Nuc_2D_unpolarized[l] - SANSData->S_NucMag_2D[l] + SANSData->S_Mag_2D_polarized[l];
		SANSData->S_Mag_2D_sanspol_p[l] = SANSData->S_Mag_2D_non_spin_flip_pp[l] + SANSData->S_Mag_2D_spin_flip_pm[l];
		SANSData->S_Mag_2D_sanspol_m[l] = SANSData->S_Mag_2D_non_spin_flip_mm[l] + SANSData->S_Mag_2D_spin_flip_mp[l];

	}

	for(unsigned long int l=0; l < (InputData->N_q); l++){

		// basic SANS cross sections
		SANSData->S_Nuc_1D_unpolarized[l] = SANSData->S_Nuc_1D_unpolarized[l] * ScalFactors->Nuc_SANS_SF;
		SANSData->S_Mag_1D_unpolarized[l] = SANSData->S_Mag_1D_unpolarized[l] * ScalFactors->Mag_SANS_SF;
		SANSData->S_Mag_1D_polarized[l] = SANSData->S_Mag_1D_polarized[l] * ScalFactors->Mag_SANS_SF;;
		SANSData->S_NucMag_1D[l] = SANSData->S_NucMag_1D[l] * ScalFactors->NucMag_SANS_SF;
		SANSData->S_Mag_1D_chiral[l] = SANSData->S_Mag_1D_chiral[l] * ScalFactors->Mag_SANS_SF;

		// composed SANS cross sections
		SANSData->S_Mag_1D_spin_flip[l] = SANSData->S_Mag_1D_unpolarized[l] - SANSData->S_Mag_1D_polarized[l];
		SANSData->S_Mag_1D_spin_flip_pm[l] = SANSData->S_Mag_1D_spin_flip[l] + SANSData->S_Mag_1D_chiral[l];
		SANSData->S_Mag_1D_spin_flip_mp[l] = SANSData->S_Mag_1D_spin_flip[l] - SANSData->S_Mag_1D_chiral[l];
		SANSData->S_Mag_1D_non_spin_flip_pp[l] = SANSData->S_Nuc_1D_unpolarized[l] + SANSData->S_NucMag_1D[l] + SANSData->S_Mag_1D_polarized[l];
		SANSData->S_Mag_1D_non_spin_flip_mm[l] = SANSData->S_Nuc_1D_unpolarized[l] - SANSData->S_NucMag_1D[l] + SANSData->S_Mag_1D_polarized[l];
		SANSData->S_Mag_1D_sanspol_p[l] = SANSData->S_Mag_1D_non_spin_flip_pp[l] + SANSData->S_Mag_1D_spin_flip_pm[l];
		SANSData->S_Mag_1D_sanspol_m[l] = SANSData->S_Mag_1D_non_spin_flip_mm[l] + SANSData->S_Mag_1D_spin_flip_mp[l];

	}


//	L = (InputData->N_r) * (InputData->N_alpha);

	for(unsigned long int l=0; l<(InputData->N_r) * (InputData->N_alpha); l++){

	    // basic Correlation Functions
	    SANSData->Corr_Nuc_2D_unpolarized[l] = SANSData->Corr_Nuc_2D_unpolarized[l] * ScalFactors->CorF_Mag_2D_SF;
		SANSData->Corr_Mag_2D_unpolarized[l] = SANSData->Corr_Mag_2D_unpolarized[l] * ScalFactors->CorF_Mag_2D_SF;
		SANSData->Corr_Mag_2D_polarized[l] = SANSData->Corr_Mag_2D_polarized[l] * ScalFactors->CorF_Mag_2D_SF;
		SANSData->Corr_NucMag_2D[l] = SANSData->Corr_NucMag_2D[l] * ScalFactors->CorF_Mag_2D_SF;
		SANSData->Corr_Mag_2D_chiral[l] = SANSData->Corr_Mag_2D_chiral[l] * ScalFactors->CorF_Mag_2D_SF;

		// composed SANS cross sections
		SANSData->Corr_Mag_2D_spin_flip[l] = SANSData->Corr_Mag_2D_unpolarized[l] - SANSData->Corr_Mag_2D_polarized[l];
		SANSData->Corr_Mag_2D_spin_flip_pm[l] = SANSData->Corr_Mag_2D_spin_flip[l] + SANSData->Corr_Mag_2D_chiral[l];
		SANSData->Corr_Mag_2D_spin_flip_mp[l] = SANSData->Corr_Mag_2D_spin_flip[l] - SANSData->Corr_Mag_2D_chiral[l];
		SANSData->Corr_Mag_2D_non_spin_flip_pp[l] = SANSData->Corr_Nuc_2D_unpolarized[l] + SANSData->Corr_NucMag_2D[l] + SANSData->Corr_Mag_2D_polarized[l];
		SANSData->Corr_Mag_2D_non_spin_flip_mm[l] = SANSData->Corr_Nuc_2D_unpolarized[l] - SANSData->Corr_NucMag_2D[l] + SANSData->Corr_Mag_2D_polarized[l];
		SANSData->Corr_Mag_2D_sanspol_p[l] = SANSData->Corr_Mag_2D_non_spin_flip_pp[l] + SANSData->Corr_Mag_2D_spin_flip_pm[l];
		SANSData->Corr_Mag_2D_sanspol_m[l] = SANSData->Corr_Mag_2D_non_spin_flip_mm[l] + SANSData->Corr_Mag_2D_spin_flip_mp[l];
		
	}

	for(unsigned long int l=0; l < (InputData->N_r); l++){

		// basic pair-distance distribution functions
		SANSData->p_Nuc_unpolarized[l] = SANSData->p_Nuc_unpolarized[l] * ScalFactors->Nuc_PairDist_SF;
		SANSData->p_Mag_unpolarized[l] = SANSData->p_Mag_unpolarized[l] * ScalFactors->Mag_PairDist_SF;
		SANSData->p_Mag_polarized[l] = SANSData->p_Mag_polarized[l] * ScalFactors->Mag_PairDist_SF;
		SANSData->p_NucMag[l] = SANSData->p_NucMag[l] * ScalFactors->NucMag_PairDist_SF;
		SANSData->p_Mag_chiral[l] = SANSData->p_Mag_chiral[l] * ScalFactors->Mag_PairDist_SF;

		// composed pair-distance distribution functions
		SANSData->p_Mag_spin_flip[l] = SANSData->p_Mag_unpolarized[l] - SANSData->p_Mag_polarized[l];
		SANSData->p_Mag_spin_flip_pm[l] = SANSData->p_Mag_spin_flip[l] + SANSData->p_Mag_chiral[l];
		SANSData->p_Mag_spin_flip_mp[l] = SANSData->p_Mag_spin_flip[l] - SANSData->p_Mag_chiral[l];
		SANSData->p_Mag_non_spin_flip_pp[l] = SANSData->p_Nuc_unpolarized[l] + SANSData->p_NucMag[l] + SANSData->p_Mag_polarized[l];
		SANSData->p_Mag_non_spin_flip_mm[l] = SANSData->p_Nuc_unpolarized[l] - SANSData->p_NucMag[l] + SANSData->p_Mag_polarized[l];
		SANSData->p_Mag_sanspol_p[l] = SANSData->p_Mag_non_spin_flip_pp[l] + SANSData->p_Mag_spin_flip_pm[l];
		SANSData->p_Mag_sanspol_m[l] = SANSData->p_Mag_non_spin_flip_mm[l] + SANSData->p_Mag_spin_flip_mp[l];



		// basic correlation functions 
		SANSData->c_Nuc_unpolarized[l] = SANSData->c_Nuc_unpolarized[l] * ScalFactors->Nuc_Corr_SF;
		SANSData->c_Mag_unpolarized[l] = SANSData->c_Mag_unpolarized[l] * ScalFactors->Mag_Corr_SF;
		SANSData->c_Mag_polarized[l] = SANSData->c_Mag_polarized[l] * ScalFactors->Mag_Corr_SF;
		SANSData->c_NucMag[l] = SANSData->c_NucMag[l] * ScalFactors->NucMag_Corr_SF;
		SANSData->c_Mag_chiral[l] = SANSData->c_Mag_chiral[l] * ScalFactors->Mag_Corr_SF;

		// composed corelation functions
		SANSData->c_Mag_spin_flip[l] = SANSData->c_Mag_unpolarized[l] - SANSData->c_Mag_polarized[l];
		SANSData->c_Mag_spin_flip_pm[l] = SANSData->c_Mag_spin_flip[l] + SANSData->c_Mag_chiral[l];
		SANSData->c_Mag_spin_flip_mp[l] = SANSData->c_Mag_spin_flip[l] - SANSData->c_Mag_chiral[l];
		SANSData->c_Mag_non_spin_flip_pp[l] = SANSData->c_Nuc_unpolarized[l] + SANSData->c_NucMag[l] + SANSData->c_Mag_polarized[l];
		SANSData->c_Mag_non_spin_flip_mm[l] = SANSData->c_Nuc_unpolarized[l] - SANSData->c_NucMag[l] + SANSData->c_Mag_polarized[l];
		SANSData->c_Mag_sanspol_p[l] = SANSData->c_Mag_non_spin_flip_pp[l] + SANSData->c_Mag_spin_flip_pm[l];
		SANSData->c_Mag_sanspol_m[l] = SANSData->c_Mag_non_spin_flip_mm[l] + SANSData->c_Mag_spin_flip_mp[l];
				
	}
		
}


// ##########################################################################################################################
void allocate_ScatteringData_RAM(InputFileData *InputData,\
						         ScatteringData *SANSData){

	SANSData->Polarization = (float*) malloc(3 * sizeof(float));
	SANSData->Polarization[0] = InputData->Polarization[0];
	SANSData->Polarization[1] = InputData->Polarization[1];
	SANSData->Polarization[2] = InputData->Polarization[2];

	SANSData->N_q = (unsigned int*) malloc(sizeof(unsigned int));
	SANSData->N_theta = (unsigned int*) malloc(sizeof(unsigned int));
	SANSData->N_r = (unsigned int*) malloc(sizeof(unsigned int));
	SANSData->N_alpha = (unsigned int*) malloc(sizeof(unsigned int));

	*SANSData->N_q = InputData->N_q;
	*SANSData->N_theta = InputData->N_theta;
	*SANSData->N_r = InputData->N_r;
	*SANSData->N_alpha = InputData->N_alpha;

	SANSData->q_max = (float*) malloc(sizeof(float));
	SANSData->r_max = (float*) malloc(sizeof(float));
	SANSData->dq = (float*) malloc(sizeof(float));
	SANSData->dtheta = (float*) malloc(sizeof(float));
	SANSData->dr = (float*) malloc(sizeof(float));
	SANSData->dalpha = (float*) malloc(sizeof(float));
	
	*SANSData->q_max = InputData->q_max;
	*SANSData->r_max = InputData->r_max;
	*SANSData->dq = (InputData->q_max)/((float) InputData->N_q - 1.0);
	*SANSData->dr = (InputData->r_max)/((float) InputData->N_r - 1.0);
	*SANSData->dtheta = (2.0*M_PI)/((float) InputData->N_theta - 1.0);
	*SANSData->dalpha = (2.0*M_PI)/((float) InputData->N_alpha - 1.0);

	unsigned int L = (InputData->N_q) * (InputData->N_theta);

	SANSData->q_2D = (float*) malloc(L*sizeof(float));
	SANSData->theta_2D = (float*) malloc(L*sizeof(float));
	SANSData->qy_2D = (float*) malloc(L*sizeof(float));
	SANSData->qz_2D = (float*) malloc(L*sizeof(float));

	SANSData->Gxx_real = (float*) malloc(L*sizeof(float));
	SANSData->Gyy_real = (float*) malloc(L*sizeof(float));
	SANSData->Gzz_real = (float*) malloc(L*sizeof(float));
	SANSData->Gxy_real = (float*) malloc(L*sizeof(float));
	SANSData->Gyx_real = (float*) malloc(L*sizeof(float));
	SANSData->Gxz_real = (float*) malloc(L*sizeof(float));
	SANSData->Gzx_real = (float*) malloc(L*sizeof(float));
	SANSData->Gyz_real = (float*) malloc(L*sizeof(float));
	SANSData->Gzy_real = (float*) malloc(L*sizeof(float));

	SANSData->Gxx_imag = (float*) malloc(L*sizeof(float));
	SANSData->Gyy_imag = (float*) malloc(L*sizeof(float));
	SANSData->Gzz_imag = (float*) malloc(L*sizeof(float));
	SANSData->Gxy_imag = (float*) malloc(L*sizeof(float));
	SANSData->Gyx_imag = (float*) malloc(L*sizeof(float));
	SANSData->Gxz_imag = (float*) malloc(L*sizeof(float));
	SANSData->Gzx_imag = (float*) malloc(L*sizeof(float));
	SANSData->Gyz_imag = (float*) malloc(L*sizeof(float));
	SANSData->Gzy_imag = (float*) malloc(L*sizeof(float));

	SANSData->S_Nuc_2D_unpolarized = (float*) malloc(L*sizeof(float));
	SANSData->S_Mag_2D_unpolarized = (float*) malloc(L*sizeof(float));
	SANSData->S_Mag_2D_polarized = (float*) malloc(L*sizeof(float));
	SANSData->S_NucMag_2D = (float*) malloc(L*sizeof(float));
	SANSData->S_Mag_2D_spin_flip = (float*) malloc(L*sizeof(float));
	SANSData->S_Mag_2D_chiral = (float*) malloc(L*sizeof(float));
	SANSData->S_Mag_2D_spin_flip_pm = (float*) malloc(L*sizeof(float));
	SANSData->S_Mag_2D_spin_flip_mp = (float*) malloc(L*sizeof(float));
	SANSData->S_Mag_2D_non_spin_flip_pp = (float*) malloc(L*sizeof(float));
	SANSData->S_Mag_2D_non_spin_flip_mm = (float*) malloc(L*sizeof(float));
	SANSData->S_Mag_2D_sanspol_p = (float*) malloc(L*sizeof(float));
	SANSData->S_Mag_2D_sanspol_m = (float*) malloc(L*sizeof(float));

	SANSData->q_1D = (float*) malloc((InputData->N_q)*sizeof(float));
	SANSData->S_Nuc_1D_unpolarized = (float*) malloc((InputData->N_q)*sizeof(float));
	SANSData->S_Mag_1D_unpolarized = (float*) malloc((InputData->N_q)*sizeof(float));
	SANSData->S_Mag_1D_polarized = (float*) malloc((InputData->N_q)*sizeof(float));
	SANSData->S_NucMag_1D = (float*) malloc((InputData->N_q)*sizeof(float));
	SANSData->S_Mag_1D_spin_flip = (float*) malloc((InputData->N_q)*sizeof(float));
	SANSData->S_Mag_1D_chiral = (float*) malloc((InputData->N_q)*sizeof(float));
	SANSData->S_Mag_1D_spin_flip_pm = (float*) malloc((InputData->N_q)*sizeof(float));
	SANSData->S_Mag_1D_spin_flip_mp = (float*) malloc((InputData->N_q)*sizeof(float));
	SANSData->S_Mag_1D_non_spin_flip_pp = (float*) malloc((InputData->N_q)*sizeof(float));
	SANSData->S_Mag_1D_non_spin_flip_mm = (float*) malloc((InputData->N_q)*sizeof(float));
	SANSData->S_Mag_1D_sanspol_p = (float*) malloc((InputData->N_q)*sizeof(float));
	SANSData->S_Mag_1D_sanspol_m = (float*) malloc((InputData->N_q)*sizeof(float));

	for(unsigned int i = 0; i < InputData->N_q; i++){
	
		SANSData->q_1D[i] = i * (*SANSData->dq);
		SANSData->S_Nuc_1D_unpolarized[i] = 0.0;
		SANSData->S_Mag_1D_unpolarized[i] = 0.0;
		SANSData->S_Mag_1D_polarized[i] = 0.0;
		SANSData->S_NucMag_1D[i] = 0.0;
		SANSData->S_Mag_1D_spin_flip[i] = 0.0;
		SANSData->S_Mag_1D_chiral[i] = 0.0;
		SANSData->S_Mag_1D_spin_flip_pm[i] = 0.0;
		SANSData->S_Mag_1D_spin_flip_mp[i] = 0.0;
		SANSData->S_Mag_1D_non_spin_flip_pp[i] = 0.0;
		SANSData->S_Mag_1D_non_spin_flip_mm[i] = 0.0;
		SANSData->S_Mag_1D_sanspol_p[i] = 0.0;
		SANSData->S_Mag_1D_sanspol_m[i] = 0.0;
		
		for(unsigned int j=0; j < InputData->N_theta; j++){
		
			SANSData->q_2D[j + i * (InputData->N_theta)] = SANSData->q_1D[i];
			SANSData->theta_2D[j + i * (InputData->N_theta)] = j * (*SANSData->dtheta);
			SANSData->qy_2D[j + i * (InputData->N_theta)] = SANSData->q_1D[i] * sin(j * (*SANSData->dtheta));
			SANSData->qz_2D[j + i * (InputData->N_theta)] = SANSData->q_1D[i] * cos(j * (*SANSData->dtheta));
			
			SANSData->Gxx_real[j + i * (InputData->N_theta)] = 0.0;
			SANSData->Gyy_real[j + i * (InputData->N_theta)] = 0.0;
			SANSData->Gzz_real[j + i * (InputData->N_theta)] = 0.0;
			SANSData->Gxy_real[j + i * (InputData->N_theta)] = 0.0;
			SANSData->Gyx_real[j + i * (InputData->N_theta)] = 0.0;
			SANSData->Gxz_real[j + i * (InputData->N_theta)] = 0.0;
			SANSData->Gzx_real[j + i * (InputData->N_theta)] = 0.0;
			SANSData->Gyz_real[j + i * (InputData->N_theta)] = 0.0;
			SANSData->Gzy_real[j + i * (InputData->N_theta)] = 0.0;

			SANSData->Gxx_imag[j + i * (InputData->N_theta)] = 0.0;
			SANSData->Gyy_imag[j + i * (InputData->N_theta)] = 0.0;
			SANSData->Gzz_imag[j + i * (InputData->N_theta)] = 0.0;
			SANSData->Gxy_imag[j + i * (InputData->N_theta)] = 0.0;
			SANSData->Gyx_imag[j + i * (InputData->N_theta)] = 0.0;
			SANSData->Gxz_imag[j + i * (InputData->N_theta)] = 0.0;
			SANSData->Gzx_imag[j + i * (InputData->N_theta)] = 0.0;
			SANSData->Gyz_imag[j + i * (InputData->N_theta)] = 0.0;
			SANSData->Gzy_imag[j + i * (InputData->N_theta)] = 0.0;

			SANSData->S_Nuc_2D_unpolarized[j + i * (InputData->N_theta)] = 0.0;
			SANSData->S_Mag_2D_unpolarized[j + i * (InputData->N_theta)] = 0.0;
			SANSData->S_NucMag_2D[j + i * (InputData->N_theta)] = 0.0;
			SANSData->S_Mag_2D_polarized[j + i * (InputData->N_theta)] = 0.0;
			SANSData->S_Mag_2D_spin_flip[j + i * (InputData->N_theta)] = 0.0;
			SANSData->S_Mag_2D_chiral[j + i * (InputData->N_theta)] = 0.0;
			SANSData->S_Mag_2D_spin_flip_pm[j + i * (InputData->N_theta)] = 0.0;
			SANSData->S_Mag_2D_spin_flip_mp[j + i * (InputData->N_theta)] = 0.0;
			SANSData->S_Mag_2D_non_spin_flip_pp[j + i * (InputData->N_theta)] = 0.0;
			SANSData->S_Mag_2D_non_spin_flip_mm[j + i * (InputData->N_theta)] = 0.0;
			SANSData->S_Mag_2D_sanspol_p[j + i * (InputData->N_theta)] = 0.0;
			SANSData->S_Mag_2D_sanspol_m[j + i * (InputData->N_theta)] = 0.0;
			
		}
	}

	L = (InputData->N_r) * (InputData->N_alpha);

	SANSData->r_2D = (float*) malloc(L*sizeof(float));
	SANSData->alpha_2D = (float*) malloc(L*sizeof(float));
	SANSData->ry_2D = (float*) malloc(L*sizeof(float));
	SANSData->rz_2D = (float*) malloc(L*sizeof(float));

	SANSData->Corr_Nuc_2D_unpolarized = (float*) malloc(L*sizeof(float));
	SANSData->Corr_Mag_2D_unpolarized = (float*) malloc(L*sizeof(float));
	SANSData->Corr_Mag_2D_polarized = (float*) malloc(L*sizeof(float));
	SANSData->Corr_NucMag_2D = (float*) malloc(L*sizeof(float));
	SANSData->Corr_Mag_2D_spin_flip = (float*) malloc(L*sizeof(float));
	SANSData->Corr_Mag_2D_chiral = (float*) malloc(L*sizeof(float));
	SANSData->Corr_Mag_2D_spin_flip_mp = (float*) malloc(L*sizeof(float));
	SANSData->Corr_Mag_2D_spin_flip_pm = (float*) malloc(L*sizeof(float));
	SANSData->Corr_Mag_2D_non_spin_flip_pp = (float*) malloc(L*sizeof(float));
	SANSData->Corr_Mag_2D_non_spin_flip_mm = (float*) malloc(L*sizeof(float));
	SANSData->Corr_Mag_2D_sanspol_p = (float*) malloc(L*sizeof(float));
	SANSData->Corr_Mag_2D_sanspol_m = (float*) malloc(L*sizeof(float));

	SANSData->r_1D = (float*) malloc((InputData->N_r)*sizeof(float));
	
	SANSData->p_Nuc_unpolarized = (float*) malloc((InputData->N_r)*sizeof(float));
	SANSData->p_Mag_unpolarized = (float*) malloc((InputData->N_r)*sizeof(float));
	SANSData->p_Mag_polarized = (float*) malloc((InputData->N_r)*sizeof(float));
	SANSData->p_NucMag = (float*) malloc((InputData->N_r)*sizeof(float));
	SANSData->p_Mag_chiral = (float*) malloc((InputData->N_r)*sizeof(float));
	SANSData->p_Mag_spin_flip = (float*) malloc((InputData->N_r)*sizeof(float));
	SANSData->p_Mag_spin_flip_pm = (float*) malloc((InputData->N_r)*sizeof(float));
	SANSData->p_Mag_spin_flip_mp = (float*) malloc((InputData->N_r)*sizeof(float));
	SANSData->p_Mag_non_spin_flip_pp = (float*) malloc((InputData->N_r)*sizeof(float));
	SANSData->p_Mag_non_spin_flip_mm = (float*) malloc((InputData->N_r)*sizeof(float));
	SANSData->p_Mag_sanspol_p = (float*) malloc((InputData->N_r)*sizeof(float));
	SANSData->p_Mag_sanspol_m = (float*) malloc((InputData->N_r)*sizeof(float));

	SANSData->c_Nuc_unpolarized = (float*) malloc((InputData->N_r)*sizeof(float));
	SANSData->c_Mag_unpolarized = (float*) malloc((InputData->N_r)*sizeof(float));
	SANSData->c_Mag_polarized = (float*) malloc((InputData->N_r)*sizeof(float));
	SANSData->c_NucMag = (float*) malloc((InputData->N_r)*sizeof(float));
	SANSData->c_Mag_chiral = (float*) malloc((InputData->N_r)*sizeof(float));
	SANSData->c_Mag_spin_flip = (float*) malloc((InputData->N_r)*sizeof(float));
	SANSData->c_Mag_spin_flip_pm = (float*) malloc((InputData->N_r)*sizeof(float));
	SANSData->c_Mag_spin_flip_mp = (float*) malloc((InputData->N_r)*sizeof(float));
	SANSData->c_Mag_non_spin_flip_pp = (float*) malloc((InputData->N_r)*sizeof(float));
	SANSData->c_Mag_non_spin_flip_mm = (float*) malloc((InputData->N_r)*sizeof(float));
	SANSData->c_Mag_sanspol_p = (float*) malloc((InputData->N_r)*sizeof(float));
	SANSData->c_Mag_sanspol_m = (float*) malloc((InputData->N_r)*sizeof(float));

	for(unsigned int i=0; i < InputData->N_r; i++){

		SANSData->r_1D[i] = i * (*SANSData->dr);
		SANSData->p_Nuc_unpolarized[i] = 0.0;
		SANSData->p_Mag_unpolarized[i] = 0.0;
		SANSData->p_Mag_polarized[i] = 0.0;
		SANSData->p_NucMag[i] = 0.0;
		SANSData->p_Mag_chiral[i] = 0.0;
		SANSData->p_Mag_spin_flip[i] = 0.0;
		SANSData->p_Mag_spin_flip_pm[i] = 0.0;
		SANSData->p_Mag_spin_flip_mp[i] = 0.0;
		SANSData->p_Mag_non_spin_flip_pp[i] = 0.0;
		SANSData->p_Mag_non_spin_flip_mm[i] = 0.0;
		SANSData->p_Mag_sanspol_p[i] = 0.0;
		SANSData->p_Mag_sanspol_m[i] = 0.0;
		
		SANSData->c_Nuc_unpolarized[i] = 0.0;
		SANSData->c_Mag_unpolarized[i] = 0.0;
		SANSData->c_Mag_polarized[i] = 0.0;
		SANSData->c_NucMag[i] = 0.0;
		SANSData->c_Mag_chiral[i] = 0.0;
		SANSData->c_Mag_spin_flip[i] = 0.0;
		SANSData->c_Mag_spin_flip_pm[i] = 0.0;
		SANSData->c_Mag_spin_flip_mp[i] = 0.0;
		SANSData->c_Mag_non_spin_flip_pp[i] = 0.0;
		SANSData->c_Mag_non_spin_flip_mm[i] = 0.0;
		SANSData->c_Mag_sanspol_p[i] = 0.0;
		SANSData->c_Mag_sanspol_m[i] = 0.0;
				
		for(unsigned int j=0; j < InputData->N_alpha; j++){

			SANSData->r_2D[j + i * (InputData->N_alpha)] = i * (*SANSData->dr);
			SANSData->alpha_2D[j + i * (InputData->N_alpha)] = j * (*SANSData->dalpha);
			SANSData->ry_2D[j + i * (InputData->N_alpha)] = SANSData->r_1D[i] * sin(j * (*SANSData->dalpha));
			SANSData->rz_2D[j + i * (InputData->N_alpha)] = SANSData->r_1D[i] * cos(j * (*SANSData->dalpha));

			SANSData->Corr_Nuc_2D_unpolarized[j + i * (InputData->N_alpha)] = 0.0;
			SANSData->Corr_Mag_2D_unpolarized[j + i * (InputData->N_alpha)] = 0.0;
			SANSData->Corr_Mag_2D_polarized[j + i * (InputData->N_alpha)] = 0.0;
			SANSData->Corr_NucMag_2D[j + i * (InputData->N_alpha)] = 0.0;
			SANSData->Corr_Mag_2D_spin_flip[j + i * (InputData->N_alpha)] = 0.0;
			SANSData->Corr_Mag_2D_chiral[j + i * (InputData->N_alpha)] = 0.0;
			SANSData->Corr_Mag_2D_spin_flip_mp[j + i * (InputData->N_alpha)] = 0.0;
			SANSData->Corr_Mag_2D_spin_flip_pm[j + i * (InputData->N_alpha)] = 0.0;
			SANSData->Corr_Mag_2D_non_spin_flip_pp[j + i * (InputData->N_alpha)] = 0.0;
			SANSData->Corr_Mag_2D_non_spin_flip_mm[j + i * (InputData->N_alpha)] = 0.0;
			SANSData->Corr_Mag_2D_sanspol_p[j + i * (InputData->N_alpha)] = 0.0;
			SANSData->Corr_Mag_2D_sanspol_m[j + i * (InputData->N_alpha)] = 0.0;
			
		}
	}
}



// ####################################################################################################################################
void allocate_ScatteringData_GPU(ScatteringData *SANSData, \
						         ScatteringData *SANSData_gpu){

	cudaMalloc(&SANSData_gpu->Polarization, 3 * sizeof(float));
	cudaMemcpy(SANSData_gpu->Polarization, SANSData->Polarization, 3*sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc(&SANSData_gpu->N_q, sizeof(unsigned int));
	cudaMalloc(&SANSData_gpu->N_theta, sizeof(unsigned int));
	cudaMalloc(&SANSData_gpu->N_r, sizeof(unsigned int));
	cudaMalloc(&SANSData_gpu->N_alpha, sizeof(unsigned int));

	cudaMemcpy(SANSData_gpu->N_q, SANSData->N_q, sizeof(unsigned int), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->N_theta, SANSData->N_theta, sizeof(unsigned int), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->N_r, SANSData->N_r, sizeof(unsigned int), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->N_alpha, SANSData->N_alpha, sizeof(unsigned int), cudaMemcpyHostToDevice);

	cudaMalloc(&SANSData_gpu->q_max, sizeof(float));
	cudaMalloc(&SANSData_gpu->r_max, sizeof(float));

	cudaMemcpy(SANSData_gpu->q_max, SANSData->q_max, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->r_max, SANSData->r_max, sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc(&SANSData_gpu->dq, sizeof(float));
	cudaMalloc(&SANSData_gpu->dtheta, sizeof(float));
	cudaMalloc(&SANSData_gpu->dr, sizeof(float));
	cudaMalloc(&SANSData_gpu->dalpha, sizeof(float));

	cudaMemcpy(SANSData_gpu->dq, SANSData->dq, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->dtheta, SANSData->dtheta, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->dr, SANSData->dr, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->dalpha, SANSData->dalpha, sizeof(float), cudaMemcpyHostToDevice);

	unsigned int L = (*SANSData->N_q) * (*SANSData->N_theta);

	cudaMalloc(&SANSData_gpu->q_2D, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->theta_2D, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->qy_2D, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->qz_2D, L*sizeof(float));

	cudaMalloc(&SANSData_gpu->Gxx_real, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Gyy_real, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Gzz_real, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Gxy_real, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Gyx_real, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Gxz_real, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Gzx_real, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Gyz_real, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Gzy_real, L*sizeof(float));

	cudaMalloc(&SANSData_gpu->Gxx_imag, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Gyy_imag, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Gzz_imag, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Gxy_imag, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Gyx_imag, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Gxz_imag, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Gzx_imag, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Gyz_imag, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Gzy_imag, L*sizeof(float));

	cudaMalloc(&SANSData_gpu->S_Nuc_2D_unpolarized, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->S_Mag_2D_unpolarized, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->S_NucMag_2D, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->S_Mag_2D_polarized, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->S_Mag_2D_spin_flip, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->S_Mag_2D_chiral, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->S_Mag_2D_spin_flip_pm, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->S_Mag_2D_spin_flip_mp, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->S_Mag_2D_non_spin_flip_pp, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->S_Mag_2D_non_spin_flip_mm, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->S_Mag_2D_sanspol_p, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->S_Mag_2D_sanspol_m, L*sizeof(float));

	cudaMalloc(&SANSData_gpu->q_1D, (*SANSData->N_q)*sizeof(float));
	cudaMalloc(&SANSData_gpu->S_Nuc_1D_unpolarized, (*SANSData->N_q)*sizeof(float));
	cudaMalloc(&SANSData_gpu->S_Mag_1D_unpolarized, (*SANSData->N_q)*sizeof(float));
	cudaMalloc(&SANSData_gpu->S_NucMag_1D, (*SANSData->N_q)*sizeof(float));
	cudaMalloc(&SANSData_gpu->S_Mag_1D_polarized, (*SANSData->N_q)*sizeof(float));
	cudaMalloc(&SANSData_gpu->S_Mag_1D_spin_flip, (*SANSData->N_q)*sizeof(float));
	cudaMalloc(&SANSData_gpu->S_Mag_1D_chiral, (*SANSData->N_q)*sizeof(float));
	cudaMalloc(&SANSData_gpu->S_Mag_1D_spin_flip_pm, (*SANSData->N_q)*sizeof(float));
	cudaMalloc(&SANSData_gpu->S_Mag_1D_spin_flip_mp, (*SANSData->N_q)*sizeof(float));
	cudaMalloc(&SANSData_gpu->S_Mag_1D_non_spin_flip_pp, (*SANSData->N_q)*sizeof(float));
	cudaMalloc(&SANSData_gpu->S_Mag_1D_non_spin_flip_mm, (*SANSData->N_q)*sizeof(float));
	cudaMalloc(&SANSData_gpu->S_Mag_1D_sanspol_p, (*SANSData->N_q)*sizeof(float));
	cudaMalloc(&SANSData_gpu->S_Mag_1D_sanspol_m, (*SANSData->N_q)*sizeof(float));

	cudaMemcpy(SANSData_gpu->q_2D, SANSData->q_2D, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->theta_2D, SANSData->theta_2D, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->qy_2D, SANSData->qy_2D, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->qz_2D, SANSData->qz_2D, L*sizeof(float), cudaMemcpyHostToDevice);

	cudaMemcpy(SANSData_gpu->Gxx_real, SANSData->Gxx_real, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Gyy_real, SANSData->Gyy_real, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Gzz_real, SANSData->Gzz_real, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Gxy_real, SANSData->Gxy_real, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Gyx_real, SANSData->Gyx_real, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Gxz_real, SANSData->Gxz_real, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Gzx_real, SANSData->Gzx_real, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Gyz_real, SANSData->Gyz_real, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Gzy_real, SANSData->Gzy_real, L*sizeof(float), cudaMemcpyHostToDevice);

	cudaMemcpy(SANSData_gpu->Gxx_imag, SANSData->Gxx_imag, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Gyy_imag, SANSData->Gyy_imag, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Gzz_imag, SANSData->Gzz_imag, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Gxy_imag, SANSData->Gxy_imag, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Gyx_imag, SANSData->Gyx_imag, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Gxz_imag, SANSData->Gxz_imag, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Gzx_imag, SANSData->Gzx_imag, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Gyz_imag, SANSData->Gyz_imag, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Gzy_imag, SANSData->Gzy_imag, L*sizeof(float), cudaMemcpyHostToDevice);

	cudaMemcpy(SANSData_gpu->S_Nuc_2D_unpolarized, SANSData->S_Nuc_2D_unpolarized, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->S_Mag_2D_unpolarized, SANSData->S_Mag_2D_unpolarized, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->S_NucMag_2D, SANSData->S_NucMag_2D, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->S_Mag_2D_polarized, SANSData->S_Mag_2D_polarized, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->S_Mag_2D_spin_flip, SANSData->S_Mag_2D_spin_flip, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->S_Mag_2D_chiral, SANSData->S_Mag_2D_chiral, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->S_Mag_2D_spin_flip_pm, SANSData->S_Mag_2D_spin_flip_pm, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->S_Mag_2D_spin_flip_mp, SANSData->S_Mag_2D_spin_flip_mp, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->S_Mag_2D_non_spin_flip_pp, SANSData->S_Mag_2D_non_spin_flip_pp, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->S_Mag_2D_non_spin_flip_mm, SANSData->S_Mag_2D_non_spin_flip_mm, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->S_Mag_2D_sanspol_p, SANSData->S_Mag_2D_sanspol_p, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->S_Mag_2D_sanspol_m, SANSData->S_Mag_2D_sanspol_m, L*sizeof(float), cudaMemcpyHostToDevice);

	cudaMemcpy(SANSData_gpu->q_1D, SANSData->q_1D, (*SANSData->N_q)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->S_Nuc_1D_unpolarized, SANSData->S_Nuc_1D_unpolarized, (*SANSData->N_q)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->S_Mag_1D_unpolarized, SANSData->S_Mag_1D_unpolarized, (*SANSData->N_q)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->S_NucMag_1D, SANSData->S_NucMag_1D, (*SANSData->N_q)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->S_Mag_1D_polarized, SANSData->S_Mag_1D_polarized, (*SANSData->N_q)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->S_Mag_1D_spin_flip, SANSData->S_Mag_1D_spin_flip, (*SANSData->N_q)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->S_Mag_1D_chiral, SANSData->S_Mag_1D_chiral, (*SANSData->N_q)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->S_Mag_1D_spin_flip_pm, SANSData->S_Mag_1D_spin_flip_pm, (*SANSData->N_q)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->S_Mag_1D_spin_flip_mp, SANSData->S_Mag_1D_spin_flip_mp, (*SANSData->N_q)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->S_Mag_1D_non_spin_flip_pp, SANSData->S_Mag_1D_non_spin_flip_pp, (*SANSData->N_q)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->S_Mag_1D_non_spin_flip_mm, SANSData->S_Mag_1D_non_spin_flip_mm, (*SANSData->N_q)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->S_Mag_1D_sanspol_p, SANSData->S_Mag_1D_sanspol_p, (*SANSData->N_q)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->S_Mag_1D_sanspol_m, SANSData->S_Mag_1D_sanspol_m, (*SANSData->N_q)*sizeof(float), cudaMemcpyHostToDevice);


	L = (*SANSData->N_r) * (*SANSData->N_alpha);

	cudaMalloc(&SANSData_gpu->r_2D, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->alpha_2D, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->ry_2D, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->rz_2D, L*sizeof(float));

	cudaMalloc(&SANSData_gpu->Corr_Nuc_2D_unpolarized, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Corr_Mag_2D_unpolarized, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Corr_Mag_2D_polarized, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Corr_NucMag_2D, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Corr_Mag_2D_spin_flip, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Corr_Mag_2D_chiral, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Corr_Mag_2D_spin_flip_mp, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Corr_Mag_2D_spin_flip_pm, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Corr_Mag_2D_non_spin_flip_pp, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Corr_Mag_2D_non_spin_flip_mm, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Corr_Mag_2D_sanspol_p, L*sizeof(float));
	cudaMalloc(&SANSData_gpu->Corr_Mag_2D_sanspol_m, L*sizeof(float));

	cudaMalloc(&SANSData_gpu->p_Nuc_unpolarized, (*SANSData->N_r)*sizeof(float));
	cudaMalloc(&SANSData_gpu->p_Mag_unpolarized, (*SANSData->N_r)*sizeof(float));
	cudaMalloc(&SANSData_gpu->p_Mag_polarized, (*SANSData->N_r)*sizeof(float));
	cudaMalloc(&SANSData_gpu->p_NucMag, (*SANSData->N_r)*sizeof(float));
	cudaMalloc(&SANSData_gpu->p_Mag_chiral, (*SANSData->N_r)*sizeof(float));
	cudaMalloc(&SANSData_gpu->p_Mag_spin_flip, (*SANSData->N_r)*sizeof(float));
	cudaMalloc(&SANSData_gpu->p_Mag_spin_flip_pm, (*SANSData->N_r)*sizeof(float));
	cudaMalloc(&SANSData_gpu->p_Mag_spin_flip_mp, (*SANSData->N_r)*sizeof(float));
	cudaMalloc(&SANSData_gpu->p_Mag_non_spin_flip_pp, (*SANSData->N_r)*sizeof(float));
	cudaMalloc(&SANSData_gpu->p_Mag_non_spin_flip_mm, (*SANSData->N_r)*sizeof(float));
	cudaMalloc(&SANSData_gpu->p_Mag_sanspol_p, (*SANSData->N_r)*sizeof(float));
	cudaMalloc(&SANSData_gpu->p_Mag_sanspol_m, (*SANSData->N_r)*sizeof(float));

	cudaMalloc(&SANSData_gpu->c_Nuc_unpolarized, (*SANSData->N_r)*sizeof(float));
	cudaMalloc(&SANSData_gpu->c_Mag_unpolarized, (*SANSData->N_r)*sizeof(float));
	cudaMalloc(&SANSData_gpu->c_Mag_polarized, (*SANSData->N_r)*sizeof(float));
	cudaMalloc(&SANSData_gpu->c_NucMag, (*SANSData->N_r)*sizeof(float));
	cudaMalloc(&SANSData_gpu->c_Mag_chiral, (*SANSData->N_r)*sizeof(float));
	cudaMalloc(&SANSData_gpu->c_Mag_spin_flip, (*SANSData->N_r)*sizeof(float));
	cudaMalloc(&SANSData_gpu->c_Mag_spin_flip_pm, (*SANSData->N_r)*sizeof(float));
	cudaMalloc(&SANSData_gpu->c_Mag_spin_flip_mp, (*SANSData->N_r)*sizeof(float));
	cudaMalloc(&SANSData_gpu->c_Mag_non_spin_flip_pp, (*SANSData->N_r)*sizeof(float));
	cudaMalloc(&SANSData_gpu->c_Mag_non_spin_flip_mm, (*SANSData->N_r)*sizeof(float));
	cudaMalloc(&SANSData_gpu->c_Mag_sanspol_p, (*SANSData->N_r)*sizeof(float));
	cudaMalloc(&SANSData_gpu->c_Mag_sanspol_m, (*SANSData->N_r)*sizeof(float));

	cudaMalloc(&SANSData_gpu->r_1D, (*SANSData->N_r)*sizeof(float));

	cudaMemcpy(SANSData_gpu->r_2D, SANSData->r_2D, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->alpha_2D, SANSData->alpha_2D, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->ry_2D, SANSData->ry_2D, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->rz_2D, SANSData->rz_2D, L*sizeof(float), cudaMemcpyHostToDevice);

	cudaMemcpy(SANSData_gpu->Corr_Nuc_2D_unpolarized, SANSData->Corr_Nuc_2D_unpolarized, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Corr_Mag_2D_unpolarized, SANSData->Corr_Mag_2D_unpolarized, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Corr_Mag_2D_polarized, SANSData->Corr_Mag_2D_polarized, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Corr_NucMag_2D, SANSData->Corr_NucMag_2D, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Corr_Mag_2D_spin_flip, SANSData->Corr_Mag_2D_spin_flip, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Corr_Mag_2D_chiral, SANSData->Corr_Mag_2D_chiral, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Corr_Mag_2D_spin_flip_mp, SANSData->Corr_Mag_2D_spin_flip_mp, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Corr_Mag_2D_spin_flip_pm, SANSData->Corr_Mag_2D_spin_flip_pm, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Corr_Mag_2D_non_spin_flip_pp, SANSData->Corr_Mag_2D_non_spin_flip_pp, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Corr_Mag_2D_non_spin_flip_mm, SANSData->Corr_Mag_2D_non_spin_flip_mm, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Corr_Mag_2D_sanspol_p, SANSData->Corr_Mag_2D_sanspol_p, L*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->Corr_Mag_2D_sanspol_m, SANSData->Corr_Mag_2D_sanspol_m, L*sizeof(float), cudaMemcpyHostToDevice);

	cudaMemcpy(SANSData_gpu->p_Nuc_unpolarized, SANSData->p_Nuc_unpolarized, (*SANSData->N_r)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->p_Mag_unpolarized, SANSData->p_Mag_unpolarized, (*SANSData->N_r)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->p_Mag_polarized, SANSData->p_Mag_polarized, (*SANSData->N_r)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->p_NucMag, SANSData->p_NucMag, (*SANSData->N_r)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->p_Mag_spin_flip, SANSData->p_Mag_spin_flip, (*SANSData->N_r)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->p_Mag_chiral, SANSData->p_Mag_chiral, (*SANSData->N_r)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->p_Mag_spin_flip_pm, SANSData->p_Mag_spin_flip_pm, (*SANSData->N_r)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->p_Mag_spin_flip_mp, SANSData->p_Mag_spin_flip_mp, (*SANSData->N_r)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->p_Mag_non_spin_flip_pp, SANSData->p_Mag_non_spin_flip_pp, (*SANSData->N_r)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->p_Mag_non_spin_flip_mm, SANSData->p_Mag_non_spin_flip_mm, (*SANSData->N_r)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->p_Mag_sanspol_p, SANSData->p_Mag_sanspol_p, (*SANSData->N_r)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->p_Mag_sanspol_m, SANSData->p_Mag_sanspol_m, (*SANSData->N_r)*sizeof(float), cudaMemcpyHostToDevice);

	cudaMemcpy(SANSData_gpu->c_Nuc_unpolarized, SANSData->c_Nuc_unpolarized, (*SANSData->N_r)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->c_Mag_unpolarized, SANSData->c_Mag_unpolarized, (*SANSData->N_r)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->c_Mag_polarized, SANSData->c_Mag_polarized, (*SANSData->N_r)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->c_NucMag, SANSData->c_NucMag, (*SANSData->N_r)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->c_Mag_spin_flip, SANSData->c_Mag_spin_flip, (*SANSData->N_r)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->c_Mag_chiral, SANSData->c_Mag_chiral, (*SANSData->N_r)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->c_Mag_spin_flip_pm, SANSData->c_Mag_spin_flip_pm, (*SANSData->N_r)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->c_Mag_spin_flip_mp, SANSData->c_Mag_spin_flip_mp, (*SANSData->N_r)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->c_Mag_non_spin_flip_pp, SANSData->c_Mag_non_spin_flip_pp, (*SANSData->N_r)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->c_Mag_non_spin_flip_mm, SANSData->c_Mag_non_spin_flip_mm, (*SANSData->N_r)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->c_Mag_sanspol_p, SANSData->c_Mag_sanspol_p, (*SANSData->N_r)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(SANSData_gpu->c_Mag_sanspol_m, SANSData->c_Mag_sanspol_m, (*SANSData->N_r)*sizeof(float), cudaMemcpyHostToDevice);

	cudaMemcpy(SANSData_gpu->r_1D, SANSData->r_1D, (*SANSData->N_r)*sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc(&SANSData_gpu, sizeof(ScatteringData));
	cudaMemcpy(SANSData_gpu, SANSData, sizeof(ScatteringData), cudaMemcpyHostToDevice);

}

// ###########################################################################################################################################
void init_ScatteringData(InputFileData *InputData,\
						 ScatteringData *SANSData, \
						 ScatteringData *SANSData_gpu){

	allocate_ScatteringData_RAM(InputData, SANSData);
	allocate_ScatteringData_GPU(SANSData, SANSData_gpu);

}

// ###########################################################################################################################################
void copyGPU2RAM_ScatteringData(ScatteringData *SANSData, \
							    ScatteringData *SANSData_gpu){


	cudaDeviceSynchronize();

	LogSystem::write("");
	LogSystem::write("copy SANSdata from GPU to RAM...");
	LogSystem::write("");

	unsigned int L = (*SANSData->N_q) * (*SANSData->N_theta);

	cudaMemcpy(SANSData->qy_2D, SANSData_gpu->qy_2D, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->qz_2D, SANSData_gpu->qz_2D, L*sizeof(float), cudaMemcpyDeviceToHost);

	cudaMemcpy(SANSData->Gxx_real, SANSData_gpu->Gxx_real, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Gyy_real, SANSData_gpu->Gyy_real, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Gzz_real, SANSData_gpu->Gzz_real, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Gxy_real, SANSData_gpu->Gxy_real, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Gyx_real, SANSData_gpu->Gyx_real, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Gxz_real, SANSData_gpu->Gxz_real, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Gzx_real, SANSData_gpu->Gzx_real, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Gyz_real, SANSData_gpu->Gyz_real, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Gzy_real, SANSData_gpu->Gzy_real, L*sizeof(float), cudaMemcpyDeviceToHost);

	cudaMemcpy(SANSData->Gxx_imag, SANSData_gpu->Gxx_imag, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Gyy_imag, SANSData_gpu->Gyy_imag, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Gzz_imag, SANSData_gpu->Gzz_imag, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Gxy_imag, SANSData_gpu->Gxy_imag, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Gyx_imag, SANSData_gpu->Gyx_imag, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Gxz_imag, SANSData_gpu->Gxz_imag, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Gzx_imag, SANSData_gpu->Gzx_imag, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Gyz_imag, SANSData_gpu->Gyz_imag, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Gzy_imag, SANSData_gpu->Gzy_imag, L*sizeof(float), cudaMemcpyDeviceToHost);

	cudaMemcpy(SANSData->S_Nuc_2D_unpolarized, SANSData_gpu->S_Nuc_2D_unpolarized, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->S_Mag_2D_unpolarized, SANSData_gpu->S_Mag_2D_unpolarized, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->S_NucMag_2D, SANSData_gpu->S_NucMag_2D, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->S_Mag_2D_polarized, SANSData_gpu->S_Mag_2D_polarized, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->S_Mag_2D_spin_flip, SANSData_gpu->S_Mag_2D_spin_flip, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->S_Mag_2D_chiral, SANSData_gpu->S_Mag_2D_chiral, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->S_Mag_2D_spin_flip_pm, SANSData_gpu->S_Mag_2D_spin_flip_pm, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->S_Mag_2D_spin_flip_mp, SANSData_gpu->S_Mag_2D_spin_flip_mp, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->S_Mag_2D_non_spin_flip_pp, SANSData_gpu->S_Mag_2D_non_spin_flip_pp, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->S_Mag_2D_non_spin_flip_mm, SANSData_gpu->S_Mag_2D_non_spin_flip_mm, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->S_Mag_2D_sanspol_p, SANSData_gpu->S_Mag_2D_sanspol_p, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->S_Mag_2D_sanspol_m, SANSData_gpu->S_Mag_2D_sanspol_m, L*sizeof(float), cudaMemcpyDeviceToHost);

	cudaMemcpy(SANSData->S_Nuc_1D_unpolarized, SANSData_gpu->S_Nuc_1D_unpolarized, (*SANSData->N_q)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->S_Mag_1D_unpolarized, SANSData_gpu->S_Mag_1D_unpolarized, (*SANSData->N_q)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->S_NucMag_1D, SANSData_gpu->S_NucMag_1D, (*SANSData->N_q)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->S_Mag_1D_polarized, SANSData_gpu->S_Mag_1D_polarized, (*SANSData->N_q)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->S_Mag_1D_spin_flip, SANSData_gpu->S_Mag_1D_spin_flip, (*SANSData->N_q)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->S_Mag_1D_chiral, SANSData_gpu->S_Mag_1D_chiral, (*SANSData->N_q)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->S_Mag_1D_spin_flip_pm, SANSData_gpu->S_Mag_1D_spin_flip_pm, (*SANSData->N_q)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->S_Mag_1D_spin_flip_mp, SANSData_gpu->S_Mag_1D_spin_flip_mp, (*SANSData->N_q)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->S_Mag_1D_non_spin_flip_pp, SANSData_gpu->S_Mag_1D_non_spin_flip_pp, (*SANSData->N_q)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->S_Mag_1D_non_spin_flip_mm, SANSData_gpu->S_Mag_1D_non_spin_flip_mm, (*SANSData->N_q)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->S_Mag_1D_sanspol_p, SANSData_gpu->S_Mag_1D_sanspol_p, (*SANSData->N_q)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->S_Mag_1D_sanspol_m, SANSData_gpu->S_Mag_1D_sanspol_m, (*SANSData->N_q)*sizeof(float), cudaMemcpyDeviceToHost);


	L = (*SANSData->N_r) * (*SANSData->N_alpha);

	cudaMemcpy(SANSData->Corr_Nuc_2D_unpolarized, SANSData_gpu->Corr_Nuc_2D_unpolarized, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Corr_Mag_2D_unpolarized, SANSData_gpu->Corr_Mag_2D_unpolarized, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Corr_Mag_2D_polarized, SANSData_gpu->Corr_Mag_2D_polarized, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Corr_NucMag_2D, SANSData_gpu->Corr_NucMag_2D, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Corr_Mag_2D_spin_flip, SANSData_gpu->Corr_Mag_2D_spin_flip, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Corr_Mag_2D_chiral, SANSData_gpu->Corr_Mag_2D_chiral, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Corr_Mag_2D_spin_flip_mp, SANSData_gpu->Corr_Mag_2D_spin_flip_mp, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Corr_Mag_2D_spin_flip_pm, SANSData_gpu->Corr_Mag_2D_spin_flip_pm, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Corr_Mag_2D_non_spin_flip_pp, SANSData_gpu->Corr_Mag_2D_non_spin_flip_pp, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Corr_Mag_2D_non_spin_flip_mm, SANSData_gpu->Corr_Mag_2D_non_spin_flip_mm, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Corr_Mag_2D_sanspol_p, SANSData_gpu->Corr_Mag_2D_sanspol_p, L*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->Corr_Mag_2D_sanspol_m, SANSData_gpu->Corr_Mag_2D_sanspol_m, L*sizeof(float), cudaMemcpyDeviceToHost);

	cudaMemcpy(SANSData->p_Nuc_unpolarized, SANSData_gpu->p_Nuc_unpolarized, (*SANSData->N_r)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->p_Mag_unpolarized, SANSData_gpu->p_Mag_unpolarized, (*SANSData->N_r)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->p_Mag_polarized, SANSData_gpu->p_Mag_polarized, (*SANSData->N_r)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->p_NucMag, SANSData_gpu->p_NucMag, (*SANSData->N_r)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->p_Mag_spin_flip, SANSData_gpu->p_Mag_spin_flip, (*SANSData->N_r)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->p_Mag_chiral, SANSData_gpu->p_Mag_chiral, (*SANSData->N_r)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->p_Mag_spin_flip_pm, SANSData_gpu->p_Mag_spin_flip_pm, (*SANSData->N_r)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->p_Mag_spin_flip_mp, SANSData_gpu->p_Mag_spin_flip_mp, (*SANSData->N_r)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->p_Mag_non_spin_flip_pp, SANSData_gpu->p_Mag_non_spin_flip_pp, (*SANSData->N_r)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->p_Mag_non_spin_flip_mm, SANSData_gpu->p_Mag_non_spin_flip_mm, (*SANSData->N_r)*sizeof(float), cudaMemcpyDeviceToHost);

	cudaMemcpy(SANSData->c_Nuc_unpolarized, SANSData_gpu->c_Nuc_unpolarized, (*SANSData->N_r)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->c_Mag_unpolarized, SANSData_gpu->c_Mag_unpolarized, (*SANSData->N_r)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->c_Mag_polarized, SANSData_gpu->c_Mag_polarized, (*SANSData->N_r)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->c_NucMag, SANSData_gpu->c_NucMag, (*SANSData->N_r)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->c_Mag_spin_flip, SANSData_gpu->c_Mag_spin_flip, (*SANSData->N_r)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->c_Mag_chiral, SANSData_gpu->c_Mag_chiral, (*SANSData->N_r)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->c_Mag_spin_flip_pm, SANSData_gpu->c_Mag_spin_flip_pm, (*SANSData->N_r)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->c_Mag_spin_flip_mp, SANSData_gpu->c_Mag_spin_flip_mp, (*SANSData->N_r)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->c_Mag_non_spin_flip_pp, SANSData_gpu->c_Mag_non_spin_flip_pp, (*SANSData->N_r)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->c_Mag_non_spin_flip_mm, SANSData_gpu->c_Mag_non_spin_flip_mm, (*SANSData->N_r)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->c_Mag_sanspol_p, SANSData_gpu->c_Mag_sanspol_p, (*SANSData->N_r)*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(SANSData->c_Mag_sanspol_m, SANSData_gpu->c_Mag_sanspol_m, (*SANSData->N_r)*sizeof(float), cudaMemcpyDeviceToHost);
	
}



void write2CSVtable_ScatteringData(InputFileData *InputData, \
					     	       ScatteringData *SANSData, \
					     	       int MagData_File_Index){

	LogSystem::write("");
	LogSystem::write("write scattering data to csv-files...");

	unsigned int L = (*SANSData->N_q) * (*SANSData->N_theta);
	
	std::string target_foldername = InputData->SANSDataFoldername + "/SANS_" + std::to_string(MagData_File_Index) + "/";
	mkdir(target_foldername.c_str(), 0777);

	ofstream fout;

	// 2D SANS cross sections

	if(InputData->output_fourier_correlation_matrix_flag || any_active(InputData->SANS2D)){
	
		std::string filename_SANS2D = target_foldername + "SANS2D.csv";	
		fout.open(filename_SANS2D);
		fout << "qz";
		fout << "," << "qy";
		fout << "," << "q";
		fout << "," << "theta";

		if(InputData->output_fourier_correlation_matrix_flag){
			fout << "," << "Gxx_real";
			fout << "," << "Gxx_imag";
		    fout << "," << "Gyy_real";
			fout << "," << "Gyy_imag";
			fout << "," << "Gzz_real"; 
			fout << "," << "Gzz_imag"; 
			fout << "," << "Gxy_real";
			fout << "," << "Gxy_imag";
			fout << "," << "Gyx_real";
			fout << "," << "Gyx_imag";
			fout << "," << "Gxz_real";
			fout << "," << "Gxz_imag";
			fout << "," << "Gzx_real";
			fout << "," << "Gzx_imag";
			fout << "," << "Gyz_real"; 
			fout << "," << "Gyz_imag";
			fout << "," << "Gzy_real"; 
			fout << "," << "Gzy_imag";
		}

		if(InputData->OutFlags.SANS2D.Nuclear){
			fout << "," << "S_N";
		}

		if(InputData->OutFlags.SANS2D.Unpolarized){
			fout << "," << "S_M";
		}

		if(InputData->OutFlags.SANS2D.NuclearMagnetic){
			fout << "," << "S_NM";
		}

		if(InputData->OutFlags.ANS2D.Polarized){
			fout << "," << "S_P";
		}

		if(InputData->OutFlags.SANS2D.Chiral){
			fout << "," << "S_chi";
		}

		if(InputData->OutFlags.SANS2D.SpinFlip){
			fout << "," << "S_sf";
		}

		if(InputData->OutFlags.SANS2D.PM_SpinFlip){
			fout << "," << "S_pm";
		}

		if(InputData->OutFlags.SANS2D.MP_SpinFlip){
			fout << "," << "S_mp";
		}

		if(InputData->OutFlags.SANS2D.PP_NonSpinFlip){
			fout << "," << "S_pp";
		}

		if(InputData->OutFlags.SANS2D.MM_NonSpinFlip){
			fout << "," << "S_mm";
		}
		
		if(InputData->OutFlags.SANS2D.P_SANSPOL){
			fout << "," << "S_p";
		}

		if(InputData->OutFlags.SANS2D.M_SANSPOL){
			fout << "," << "S_m";
		}

		fout << "\n";
		
		for(unsigned long int n = 0; n<L; n++){
			fout << SANSData->qz_2D[n];
			fout << "," << SANSData->qy_2D[n];
			fout << "," << SANSData->q_2D[n];
			fout << "," << SANSData->theta_2D[n];
			
			if(InputData->output_fourier_correlation_matrix_flag){
				fout << "," << SANSData->Gxx_real[n];
				fout << "," << SANSData->Gxx_imag[n];
				fout << "," << SANSData->Gyy_real[n];
				fout << "," << SANSData->Gyy_imag[n];
				fout << "," << SANSData->Gzz_real[n];
				fout << "," << SANSData->Gzz_imag[n];
				fout << "," << SANSData->Gxy_real[n];
				fout << "," << SANSData->Gxy_imag[n];
				fout << "," << SANSData->Gyx_real[n];
				fout << "," << SANSData->Gyx_imag[n];
				fout << "," << SANSData->Gxz_real[n];
				fout << "," << SANSData->Gxz_imag[n];
				fout << "," << SANSData->Gzx_real[n];
				fout << "," << SANSData->Gzx_imag[n];
				fout << "," << SANSData->Gyz_real[n];
				fout << "," << SANSData->Gyz_imag[n];
				fout << "," << SANSData->Gzy_real[n];
				fout << "," << SANSData->Gzy_imag[n];
			}


			if(InputData->OutFlags.SANS2D.Nuclear){
				fout << "," << SANSData->S_Nuc_2D_unpolarized[n];
			}

			if(InputData->OutFlags.SANS2D.Unpolarized){
				fout << "," << SANSData->S_Mag_2D_unpolarized[n];
			}

			if(InputData->OutFlags.SANS2D.NuclearMagnetic){
				fout << "," << SANSData->S_NucMag_2D[n];
			}
			
			if(InputData->OutFlags.SANS2D.Polarized){
				fout << "," << SANSData->S_Mag_2D_polarized[n];
			}

			if(InputData->OutFlags.SANS2D.Chiral){
				fout << "," << SANSData->S_Mag_2D_chiral[n];
			}
			
			if(InputData->OutFlags.SANS2D.SpinFlip){
				fout << "," << SANSData->S_Mag_2D_spin_flip[n];
			}

			if(InputData->OutFlags.SANS2D.PM_SpinFlip){
				fout << "," << SANSData->S_Mag_2D_spin_flip_pm[n];
			}
		
			if(InputData->OutFlags.SANS2D.MP_SpinFlip){
				fout << "," << SANSData->S_Mag_2D_spin_flip_mp[n];
			}
		
			if(InputData->OutFlags.SANS2D.PP_NonSpinFlip){
				fout << "," << SANSData->S_Mag_2D_non_spin_flip_pp[n];
			}
		
			if(InputData->OutFlags.SANS2D.MM_NonSpinFlip){
				fout << "," << SANSData->S_Mag_2D_non_spin_flip_mm[n];
			}
			
			if(InputData->OutFlags.SANS2D.P_SANSPOL){
				fout << "," << SANSData->S_Mag_2D_sanspol_p[n];
			}

			if(InputData->OutFlags.SANS2D.M_SANSPOL){
				fout << "," << SANSData->S_Mag_2D_sanspol_m[n];
			}
			

			fout << "\n";
		}
		fout.close();
		LogSystem::write("SANS2D finished...");
	}

	// 1D SANS cross sections ####################################################
	if(any_active(InputData->OutFlags.SANS1D)){
	
	std::string filename_SANS1D = target_foldername + "SANS1D.csv";
	fout.open(filename_SANS1D);

	fout << "q"; 

	if(InputData->OutFlags.SANS1D.Nuclear){
		fout << "," << "I_N";
	}

	if(InputData->OutFlags.SANS1D.Unpolarized){
		fout << "," << "I_M";
	}

	if(InputData->OutFlags.SANS1D.NuclearMagnetic){
		fout << "," << "I_NM";
	}

	if(InputData->OutFlags.SANS1D.Polarized){
		fout << "," << "I_P";
	}

	if(InputData->OutFlags.SANS1D.Chiral){
		fout << "," << "I_chi";
	}

	if(InputData->OutFlags.SANS1D.SpinFlip){
		fout << "," << "I_sf";
	}

	if(InputData->OutFlags.SANS1D.PM_SpinFlip){
		fout << "," << "I_pm";
	}

	if(InputData->OutFlags.SANS1D.MP_SpinFlip){
		fout << "," << "I_mp";
	}

	if(InputData->OutFlags.SANS1D.PP_NonSpinFlip){
		fout << "," << "I_pp";
	}

	if(InputData->OutFlags.SANS1D.MM_NonSpinFlip){
		fout << "," << "I_mm";
	}

	if(InputData->OutFlags.SANS1D.P_SANSPOL){
		fout << "," << "I_p";
	}

	if(InputData->OutFlags.SANS1D.M_SANSPOL){
		fout << "," << "I_m";
	}

	fout << "\n";	

	for(unsigned long int n = 0; n<(*SANSData->N_q); n++){

		fout << SANSData->q_1D[n];

		if(InputData->OutFlags.SANS1D.Nuclear){
			fout << "," << SANSData->S_Nuc_1D_unpolarized[n];
		}
	
		if(InputData->OutFlags.SANS1D.Unpolarized){
			fout << "," << SANSData->S_Mag_1D_unpolarized[n];
		}

		if(InputData->OutFlags.SANS1D.NuclearMagnetic){
			fout << "," << SANSData->S_NucMag_1D[n];
		}
	
		if(InputData->OutFlags.SANS1D.Polarized){
			fout << "," << SANSData->S_Mag_1D_polarized[n];
		}
	
		if(InputData->OutFlags.SANS1D.Chiral){
			fout << "," << SANSData->S_Mag_1D_chiral[n];
		}	
	
		if(InputData->OutFlags.SANS1D.SpinFlip){
			fout << "," << SANSData->S_Mag_1D_spin_flip[n];
		}
	
		if(InputData->OutFlags.SANS1D.PM_SpinFlip){
			fout << "," << SANSData->S_Mag_1D_spin_flip_pm[n];
		}

		if(InputData->OutFlags.SANS1D.MP_SpinFlip){
			fout << "," << SANSData->S_Mag_1D_spin_flip_mp[n];
		}

		if(InputData->OutFlags.SANS1D.PP_NonSpinFlip){
			fout << "," << SANSData->S_Mag_1D_non_spin_flip_pp[n];
		}

		if(InputData->OutFlags.SANS1D.MM_NonSpinFlip){
			fout << "," << SANSData->S_Mag_1D_non_spin_flip_mm[n];
		}

		if(InputData->OutFlags.SANS1D.P_SANSPOL){
			fout << "," << SANSData->S_Mag_1D_sanspol_p[n];
		}

		if(InputData->OutFlags.SANS1D.M_SANSPOL){
			fout << "," << SANSData->S_Mag_1D_sanspol_m[n];
		}

		fout << "\n";
	}
	fout.close();
	LogSystem::write("SANS1D finished...");
	}
	

	// Pair-Distance Distribution 1D ######################################################################

	if(any_active(InputData->OutFlags.Corr1D) || any_active(InputData->OutFlags.PairDist1D)){

		std::string filename_Corr1D = target_foldername + "Corr1D.csv";
		fout.open(filename_Corr1D);

		fout << "r";

		if(InputData->OutFlags.PairDist1D.Nuclear){
			fout << "," << "p_N";
		}

		if(InputData->OutFlags.PairDist1D.Unpolarized){
			fout << "," << "p_M";
		}

		if(InputData->OutFlags.PairDist1D.NuclearMagnetic){
			fout << "," << "p_NM";
		}

		if(InputData->OutFlags.PairDist1D.Polarized){
			fout << "," << "p_P";
		}

		if(InputData->OutFlags.PairDist1D.Chiral){
			fout << "," << "p_chi";
		}

		if(InputData->OutFlags.PairDist1D.SpinFlip){
			fout << "," << "p_sf";
		}

		if(InputData->OutFlags.PairDist1D.PM_SpinFlip){
			fout << "," << "p_pm";
		}

		if(InputData->OutFlags.PairDist1D.MP_SpinFlip){
			fout << "," << "p_mp";
		}

		if(InputData->OutFlags.PairDist1D.PP_NonSpinFlip){
			fout << "," << "p_pp";
		}

		if(InputData->OutFlags.PairDist1D.MM_NonSpinFlip){
			fout << "," << "p_mm";
		}
		
		if(InputData->OutFlags.PairDist1D.P_SANSPOL){
			fout << "," << "p_p";
		}

		if(InputData->OutFlags.PairDist1D.M_SANSPOL){
			fout << "," << "p_m";
		}


		

		if(InputData->OutFlags.Corr1D.Nuclear){
			fout << "," << "c_N";
		}

		if(InputData->OutFlags.Corr1D.Unpolarized){
			fout << "," << "c_M";
		}

		if(InputData->OutFlags.Corr1D.NuclearMagnetic){
			fout << "," << "c_NM";
		}

		if(InputData->OutFlags.Corr1D.Polarized){
			fout << "," << "c_P";
		}

		if(InputData->OutFlags.Corr1D.Chiral){
			fout << "," << "c_chi";
		}

		if(InputData->OutFlags.Corr1D.SpinFlip){
			fout << "," << "c_sf";
		}

		if(InputData->OutFlags.Corr1D.PM_SpinFlip){
			fout << "," << "c_pm";
		}

		if(InputData->OutFlags.Corr1D.MP_SpinFlip){
			fout << "," << "c_mp";
		}

		if(InputData->OutFlags.Corr1D.PP_SpinFlip){
			fout << "," << "c_pp";
		}

		if(InputData->OutFlags.Corr1D.MM_NonSpinFlip){
			fout << "," << "c_mm";
		}
		
		if(InputData->OutFlags.Corr1D.P_SANSPOL){
			fout << "," << "c_p";
		}

		if(InputData->OutFlags.Corr1D.M_SANSPOL){
			fout << "," << "c_m";
		}
		
		fout << "\n";

		for(unsigned long int n = 0; n<(*SANSData->N_r); n++){
		
			fout << SANSData->r_1D[n];

			if(InputData->OutFlags.PairDist1D.Nuclear){
				fout << "," << SANSData->p_Nuc_unpolarized[n];
			}

			if(InputData->OutFlags.PairDist1D.Unpolarized){
				fout << "," << SANSData->p_Mag_unpolarized[n];
			}

			if(InputData->OutFlags.PairDist1D.NuclearMagnetic){
				fout << "," << SANSData->p_NucMag[n];
			}

			if(InputData->OutFlags.PairDist1D.Polarized){
				fout << "," << SANSData->p_Mag_polarized[n];
			}

			if(InputData->OutFlags.PairDist1D.Chiral){
				fout << "," << SANSData->p_Mag_chiral[n];
			}

			if(InputData->OutFlags.PairDist1D.SpinFlip){
				fout << "," << SANSData->p_Mag_spin_flip[n];
			}

			if(InputData->OutFlags.PairDist1D.PM_SpinFlip){
				fout << "," << SANSData->p_Mag_spin_flip_pm[n];
			}
	
			if(InputData->OutFlags.PairDist1D.MP_SpinFlip){
				fout << "," << SANSData->p_Mag_spin_flip_mp[n];
			}
	
			if(InputData->OutFlags.PairDist1D.PP_NonSpinFlip){
				fout << "," << SANSData->p_Mag_non_spin_flip_pp[n];
			}
	
			if(InputData->OutFlags.PairDist1D.MM_NonSpinFlip){
				fout << "," << SANSData->p_Mag_non_spin_flip_mm[n];
			}
			
			if(InputData->OutFlags.PairDist1D.P_SANSPOL){
				fout << "," << SANSData->p_Mag_sanspol_p[n];
			}
	
			if(InputData->OutFlags.PairDist1D.M_SANSPOL){
				fout << "," << SANSData->p_Mag_sanspol_m[n];
			}


			

			if(InputData->OutFlags.Corr1D.Nuclear){
				fout << "," << SANSData->c_Nuc_unpolarized[n];
			}

			if(InputData->OutFlags.Corr1D.Unpolarized){
				fout << "," << SANSData->c_Mag_unpolarized[n];
			}

			if(InputData->OutFlags.Corr1D.NuclearMagnetic){
				fout << "," << SANSData->c_NucMag[n];
			}

			if(InputData->OutFlags.Corr1D.Polarized){
				fout << "," << SANSData->c_Mag_polarized[n];
			}

			if(InputData->OutFlags.Corr1D.Chiral){
				fout << "," << SANSData->c_Mag_chiral[n];
			}

			if(InputData->OutFlags.Corr1D.SpinFlip){
				fout << "," << SANSData->c_Mag_spin_flip[n];
			}

			if(InputData->OutFlags.Corr1D.PM_SpinFlip){
				fout << "," << SANSData->c_Mag_spin_flip_pm[n];
			}
	
			if(InputData->OutFlags.Corr1D.MP_SpinFlip){
				fout << "," << SANSData->c_Mag_spin_flip_mp[n];
			}
	
			if(InputData->OutFlags.Corr1D.PP_NonSpinFlip){
				fout << "," << SANSData->c_Mag_non_spin_flip_pp[n];
			}
	
			if(InputData->OutFlags.Corr1D.MM_NonSpinFlip){
				fout << "," << SANSData->c_Mag_non_spin_flip_mm[n];
			}
			
			if(InputData->OutFlags.Corr1D.P_SANSPOL){
				fout << "," << SANSData->c_Mag_sanspol_p[n];
			}
	
			if(InputData->OutFlags.Corr1D.M_SANSPOL){
				fout << "," << SANSData->c_Mag_sanspol_m[n];
			}

			fout << "\n";
			
		}
		fout.close();
		LogSystem::write("Corr1D finished...");
	}

	// Corr 2D ###################################################################################################################################

	if(any_active(InputData->OutFlags.Corr2D)){

		L = (*SANSData->N_r) * (*SANSData->N_alpha);

		std::string filename_Corr2D = target_foldername + "Corr2D.csv";
		ofstream fout;
		fout.open(filename_Corr2D);
		fout << "rz";
		fout << "," << "ry";
		fout << "," << "r";
		fout << "," << "alpha";

		if(InputData->OutFlags.Corr2D.Nuclear){
			fout << "," << "C_N";
		}

		if(InputData->OutFlags.Corr2D.Unpolarized){
			fout << "," << "C_M";
		}

		if(InputData->OutFlags.Corr2D.NuclearMagnetic){
			fout << "," << "C_NM";
		}

		if(InputData->OutFlags.Corr2D.Polarized){
			fout << "," << "C_P";
		}

		if(InputData->OutFlags.Corr2D.Chiral){
			fout << "," << "C_chi";
		}

		if(InputData->OutFlags.Corr2D.SpinFlip){
			fout << "," << "C_sf";
		}

		if(InputData->OutFlags.Corr2D.PM_SpinFlip){
			fout << "," << "C_pm";
		}

		if(InputData->OutFlags.Corr2D.MP_SpinFlip){
			fout << "," << "C_mp";
		}

		if(InputData->OutFlags.Corr2D.PP_NonSpinFlip){
			fout << "," << "C_pp";
		}

		if(InputData->OutFlags.Corr2D.MM_NonSpinFlip){
			fout << "," << "C_mm";
		}

		if(InputData->OutFlags.Corr2D.P_SANSPOL){
			fout << "," << "C_p";
		}

		if(InputData->OutFlags.Corr2D.M_SANSPOL){
			fout << "," << "C_m";
		}

		fout << "\n";
		
		for(unsigned long int n=0; n<L; n++){

			fout << SANSData->rz_2D[n];
			fout << "," << SANSData->ry_2D[n];
			fout << "," << SANSData->r_2D[n];
			fout << "," << SANSData->alpha_2D[n];

			if(InputData->OutFlags.Corr2D.Nuclear){
				fout << "," << SANSData->Corr_Nuc_2D_unpolarized[n];
			}

			if(InputData->OutFlags.Corr2D.Unpolarized){
				fout << "," << SANSData->Corr_Mag_2D_unpolarized[n];
			}

			if(InputData->OutFlags.Corr2D.NuclearMagnetic){
				fout << "," << SANSData->Corr_NucMag_2D[n];
			}

			if(InputData->OutFlags.Corr2D.Polarized){
				fout << "," << SANSData->Corr_Mag_2D_polarized[n];
			}

			if(InputData->OutFlags.Corr2D.Chiral){
				fout << "," << SANSData->Corr_Mag_2D_chiral[n];
			}
		
			if(InputData->OutFlags.Corr2D.SpinFlip){
				fout << "," << SANSData->Corr_Mag_2D_spin_flip[n];
			}

			if(InputData->OutFlags.Corr2D.PM_SpinFlip){
				fout << "," << SANSData->Corr_Mag_2D_spin_flip_pm[n];
			}

			if(InputData->OutFlags.Corr2D.MP_SpinFlip){
				fout << "," << SANSData->Corr_Mag_2D_spin_flip_mp[n];
			}

			if(InputData->OutFlags.Corr2D.PP_NonSpinFlip){
				fout << "," << SANSData->Corr_Mag_2D_non_spin_flip_pp[n];
			}

			if(InputData->OutFlags.Corr2D.MM_NonSpinFlip){
				fout << "," << SANSData->Corr_Mag_2D_non_spin_flip_mm[n];
			}

			if(InputData->OutFlags.Corr2D.P_SANSPOL){
				fout << "," << SANSData->Corr_Mag_2D_sanspol_p[n];
			}

			if(InputData->OutFlags.Corr2D.M_SANSPOL){
				fout << "," << SANSData->Corr_Mag_2D_sanspol_m[n];
			}

			fout << "\n";
			
		}
		fout.close();
		LogSystem::write("Corr2D finished...");
	}
	
}



// struct Column2D {
//     const char* header;
//     const double* data;
//     bool enabled;
// };

// void write2CSVtable_ScatteringData(InputFileData *InputData, \
//  					     	       ScatteringData *SANSData, \
// 					     	       int MagData_File_Index)
// {

// 	std::vector<Column2D> columns;

// 	if (InputData->OutFlags.SANS2D.Nuclear)
//     	columns.push_back({"S_N", SANSData->S_Nuc_2D_unpolarized, true});

// 	if (InputData->OutFlags.SANS2D.Unpolarized)
// 		columns.push_back({"S_M", SANSData->S_Mag_2D_unpolarized, true});

// 	if (InputData->OutFlags.SANS2D.NuclearMagnetic)
// 		columns.push_back({"S_NM", SANSData->S_NucMag_2D, true});

// 	if (InputData->OutFlags.SANS2D.SpinFlip)
// 		columns.push_back({"S_sf", SANSData->S_Mag_2D_spin_flip, true});

// 	if (InputData->OutFlags.SANS2D.Chiral)
// 		columns.push_back({"S_chi", SANSData->S_Mag_2D_chiral, true});

// 	if (InputData->OutFlags.SANS2D.PM_SpinFlip)
// 		columns.push_back({"S_pm", SANSData->S_Mag_2D_spin_flip_pm, true});

// 	if (InputData->OutFlags.SANS2D.MP_SpinFlip)
// 		columns.push_back({"S_mp", SANSData->S_Mag_2D_spin_flip_mp, true});

// 	if (InputData->OutFlags.SANS2D.PP_NonSpinFlip)
// 		columns.push_back({"S_pp", SANSData->S_Mag_2D_non_spin_flip_pp, true});

// 	if (InputData->OutFlags.SANS2D.MM_NonSpinFlip)
// 		columns.push_back({"S_mm", SANSData->S_Mag_2D_non_spin_flip_mm, true});

// 	if (InputData->OutFlags.SANS2D.P_SANSPOL)
// 		columns.push_back({"S_mm", SANSData->S_Mag_2D_sanspol_p, true});

// 	if (InputData->OutFlags.SANS2D.M_SANSPOL)
// 		columns.push_back({"S_mm", SANSData->S_Mag_2D_sanspol_M, true});


// 	fout << "qz,qy,q,theta";

// 	for (const auto& col : columns)
// 		fout << "," << col.header;

// 	fout << "\n";

// 	for (size_t n = 0; n < L; ++n) {

//     fout << SANSData->qz_2D[n]
//          << "," << SANSData->qy_2D[n]
//          << "," << SANSData->q_2D[n]
//          << "," << SANSData->theta_2D[n];

//     for (const auto& col : columns)
//         fout << "," << col.data[n];

//     fout << "\n";
// }



// }


void free_ScatteringData(ScatteringData *SANSData, \
				         ScatteringData *SANSData_gpu){

	LogSystem::write("free SANSdata...");

//	cudaError_t err;

	free(SANSData->Polarization);

	free(SANSData->N_q);
	free(SANSData->N_theta);
	free(SANSData->N_r);
	free(SANSData->N_alpha);

	free(SANSData->q_max);
	free(SANSData->r_max);
	free(SANSData->dq);
	free(SANSData->dtheta);
	free(SANSData->dr);
	free(SANSData->dalpha);

	free(SANSData->qy_2D);
	free(SANSData->qz_2D);
	free(SANSData->q_2D);
	free(SANSData->theta_2D);

	free(SANSData->ry_2D);
	free(SANSData->rz_2D);
	free(SANSData->r_2D);
	free(SANSData->alpha_2D);

	free(SANSData->q_1D);
	free(SANSData->r_1D);

	free(SANSData->Gxx_real);
	free(SANSData->Gyy_real);
	free(SANSData->Gzz_real);
	free(SANSData->Gxy_real);
	free(SANSData->Gyx_real);
	free(SANSData->Gxz_real);
	free(SANSData->Gzx_real);
	free(SANSData->Gyz_real);
	free(SANSData->Gzy_real);

	free(SANSData->Gxx_imag);
	free(SANSData->Gyy_imag);
	free(SANSData->Gzz_imag);
	free(SANSData->Gxy_imag);
	free(SANSData->Gyx_imag);
	free(SANSData->Gxz_imag);
	free(SANSData->Gzx_imag);
	free(SANSData->Gyz_imag);
	free(SANSData->Gzy_imag);

	free(SANSData->S_Nuc_2D_unpolarized);
	free(SANSData->S_Mag_2D_unpolarized);
	free(SANSData->S_NucMag_2D);
	free(SANSData->S_Mag_2D_polarized);
	free(SANSData->S_Mag_2D_spin_flip);
	free(SANSData->S_Mag_2D_chiral);
	free(SANSData->S_Mag_2D_spin_flip_pm);
	free(SANSData->S_Mag_2D_spin_flip_mp);
	free(SANSData->S_Mag_2D_non_spin_flip_pp);
	free(SANSData->S_Mag_2D_non_spin_flip_mm);
	free(SANSData->S_Mag_2D_sanspol_p);
	free(SANSData->S_Mag_2D_sanspol_m);

	free(SANSData->Corr_Nuc_2D_unpolarized);
	free(SANSData->Corr_Mag_2D_unpolarized);
	free(SANSData->Corr_Mag_2D_polarized);
	free(SANSData->Corr_NucMag_2D);
	free(SANSData->Corr_Mag_2D_spin_flip);
	free(SANSData->Corr_Mag_2D_chiral);
	free(SANSData->Corr_Mag_2D_spin_flip_mp);
	free(SANSData->Corr_Mag_2D_spin_flip_pm);
	free(SANSData->Corr_Mag_2D_non_spin_flip_pp);
	free(SANSData->Corr_Mag_2D_non_spin_flip_mm);
	free(SANSData->Corr_Mag_2D_sanspol_p);
	free(SANSData->Corr_Mag_2D_sanspol_m);

	free(SANSData->S_Nuc_1D_unpolarized);
	free(SANSData->S_Mag_1D_unpolarized);
	free(SANSData->S_NucMag_1D);
	free(SANSData->S_Mag_1D_polarized);
	free(SANSData->S_Mag_1D_chiral);
	free(SANSData->S_Mag_1D_spin_flip);
	free(SANSData->S_Mag_1D_spin_flip_pm);
	free(SANSData->S_Mag_1D_spin_flip_mp);
	free(SANSData->S_Mag_1D_non_spin_flip_pp);
	free(SANSData->S_Mag_1D_non_spin_flip_mm);
	free(SANSData->S_Mag_1D_sanspol_p);
	free(SANSData->S_Mag_1D_sanspol_m);

	free(SANSData->p_Nuc_unpolarized);
	free(SANSData->p_Mag_unpolarized);
	free(SANSData->p_Mag_polarized);
	free(SANSData->p_NucMag);
	free(SANSData->p_Mag_chiral);
	free(SANSData->p_Mag_spin_flip);
	free(SANSData->p_Mag_spin_flip_pm);
	free(SANSData->p_Mag_spin_flip_mp);
	free(SANSData->p_Mag_non_spin_flip_pp);
	free(SANSData->p_Mag_non_spin_flip_mm);
	free(SANSData->p_Mag_sanspol_p);
	free(SANSData->p_Mag_sanspol_m);

	free(SANSData->c_Nuc_unpolarized);
	free(SANSData->c_Mag_unpolarized);
	free(SANSData->c_Mag_polarized);
	free(SANSData->c_NucMag);
	free(SANSData->c_Mag_chiral);
	free(SANSData->c_Mag_spin_flip);
	free(SANSData->c_Mag_spin_flip_pm);
	free(SANSData->c_Mag_spin_flip_mp);
	free(SANSData->c_Mag_non_spin_flip_pp);
	free(SANSData->c_Mag_non_spin_flip_mm);
	free(SANSData->c_Mag_sanspol_p);
	free(SANSData->c_Mag_sanspol_m);

	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Polarization);
	cudaFree(SANSData_gpu->N_q);
	cudaFree(SANSData_gpu->N_theta);
	cudaFree(SANSData_gpu->N_r);
	cudaFree(SANSData_gpu->N_alpha);
	cudaDeviceSynchronize();

	cudaFree(SANSData_gpu->q_max);
	cudaFree(SANSData_gpu->r_max);
	cudaFree(SANSData_gpu->dq);
	cudaFree(SANSData_gpu->dtheta);
	cudaFree(SANSData_gpu->dr);
	cudaFree(SANSData_gpu->dalpha);
	cudaDeviceSynchronize();

	cudaFree(SANSData_gpu->qy_2D);
	cudaFree(SANSData_gpu->qz_2D);
	cudaFree(SANSData_gpu->q_2D);
	cudaFree(SANSData_gpu->theta_2D);
	cudaDeviceSynchronize();

	cudaFree(SANSData_gpu->ry_2D);
	cudaFree(SANSData_gpu->rz_2D);
	cudaFree(SANSData_gpu->r_2D);
	cudaFree(SANSData_gpu->alpha_2D);
	cudaDeviceSynchronize();

	cudaFree(SANSData_gpu->q_1D);
	cudaFree(SANSData_gpu->r_1D);
	cudaDeviceSynchronize();

	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Gxx_real);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Gyy_real);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Gzz_real);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Gxy_real);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Gyx_real);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Gxz_real);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Gzx_real);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Gyz_real);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Gzy_real);
	cudaDeviceSynchronize();

	cudaFree(SANSData_gpu->Gxx_imag);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Gyy_imag);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Gzz_imag);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Gxy_imag);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Gyx_imag);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Gxz_imag);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Gzx_imag);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Gyz_imag);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Gzy_imag);
	cudaDeviceSynchronize();

	cudaFree(SANSData_gpu->S_Nuc_2D_unpolarized);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->S_Mag_2D_unpolarized);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->S_NucMag_2D);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->S_Mag_2D_polarized);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->S_Mag_2D_spin_flip);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->S_Mag_2D_chiral);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->S_Mag_2D_spin_flip_pm);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->S_Mag_2D_spin_flip_mp);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->S_Mag_2D_non_spin_flip_pp);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->S_Mag_2D_non_spin_flip_mm);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->S_Mag_2D_sanspol_p);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->S_Mag_2D_sanspol_m);
	cudaDeviceSynchronize();

	cudaFree(SANSData_gpu->Corr_Nuc_2D_unpolarized);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Corr_Mag_2D_unpolarized);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Corr_Mag_2D_polarized);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Corr_NucMag_2D);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Corr_Mag_2D_spin_flip);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Corr_Mag_2D_chiral);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Corr_Mag_2D_spin_flip_mp);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Corr_Mag_2D_spin_flip_pm);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Corr_Mag_2D_non_spin_flip_pp);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Corr_Mag_2D_non_spin_flip_mm);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Corr_Mag_2D_sanspol_p);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->Corr_Mag_2D_sanspol_m);
	cudaDeviceSynchronize();

	cudaFree(SANSData_gpu->S_Nuc_1D_unpolarized);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->S_Mag_1D_unpolarized);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->S_NucMag_1D);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->S_Mag_1D_polarized);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->S_Mag_1D_chiral);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->S_Mag_1D_spin_flip);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->S_Mag_1D_spin_flip_pm);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->S_Mag_1D_spin_flip_mp);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->S_Mag_1D_non_spin_flip_pp);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->S_Mag_1D_non_spin_flip_mm);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->S_Mag_1D_sanspol_p);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->S_Mag_1D_sanspol_m);
	cudaDeviceSynchronize();


	cudaFree(SANSData_gpu->p_Nuc_unpolarized);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->p_Mag_unpolarized);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->p_Mag_polarized);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->p_NucMag);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->p_Mag_chiral);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->p_Mag_spin_flip);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->p_Mag_spin_flip_pm);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->p_Mag_spin_flip_mp);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->p_Mag_non_spin_flip_pp);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->p_Mag_non_spin_flip_mm);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->p_Mag_sanspol_p);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->p_Mag_sanspol_m);
	cudaDeviceSynchronize();

	cudaFree(SANSData_gpu->c_Nuc_unpolarized);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->c_Mag_unpolarized);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->c_Mag_polarized);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->c_NucMag);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->c_Mag_chiral);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->c_Mag_spin_flip);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->c_Mag_spin_flip_pm);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->c_Mag_spin_flip_mp);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->c_Mag_non_spin_flip_pp);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->c_Mag_non_spin_flip_mm);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->c_Mag_sanspol_p);
	cudaDeviceSynchronize();
	cudaFree(SANSData_gpu->c_Mag_sanspol_m);
	cudaDeviceSynchronize();

	LogSystem::write("free SANSdata finished.");

}
