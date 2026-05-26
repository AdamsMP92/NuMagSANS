#pragma once

// Computes combined nuclear and magnetic atomistic SANS for a dilute
// orientation ensemble.
// Active data layers:
//   - NuclearData: local nuclear scattering length density data.
//   - MagnetizationData: local magnetic moments.
//   - RotationData: object-wise local orientation rotations.
// Inactive data layers:
//   - StructureData
// Use case:
//   Nuclear-magnetic scattering from a dilute ensemble of rotated local objects
//   without explicit object-center positions.
//
// Placeholder for:
// Atomistic_NuMagSANS_Kernel_RotDilute


__global__
void Atomistic_NuMagSANS_Kernel_RotDilute(NuclearData NucData,\
									   MagnetizationData MagData,\
                                       RotationData RotData,\
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

	float xr = 0.0;
	float yr = 0.0;
	float zr = 0.0;
	float mxr = 0.0;
	float myr = 0.0;
	float mzr = 0.0;
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

				// individual position rotation
				xr = RotData.RotMat[9*k+0] * MagData.x[l+N_cum] \
				   + RotData.RotMat[9*k+3] * MagData.y[l+N_cum] \
				   + RotData.RotMat[9*k+6] * MagData.z[l+N_cum];
				yr = RotData.RotMat[9*k+1] * MagData.x[l+N_cum] \
				   + RotData.RotMat[9*k+4] * MagData.y[l+N_cum] \
				   + RotData.RotMat[9*k+7] * MagData.z[l+N_cum];
				zr = RotData.RotMat[9*k+2] * MagData.x[l+N_cum] \
				   + RotData.RotMat[9*k+5] * MagData.y[l+N_cum] \
				   + RotData.RotMat[9*k+8] * MagData.z[l+N_cum];

				mxr = RotData.RotMat[9*k+0] * MagData.mx[l+N_cum] \
				   + RotData.RotMat[9*k+3] * MagData.my[l+N_cum] \
				   + RotData.RotMat[9*k+6] * MagData.mz[l+N_cum];
				myr = RotData.RotMat[9*k+1] * MagData.mx[l+N_cum] \
				   + RotData.RotMat[9*k+4] * MagData.my[l+N_cum] \
				   + RotData.RotMat[9*k+7] * MagData.mz[l+N_cum];
				mzr = RotData.RotMat[9*k+2] * MagData.mx[l+N_cum] \
				   + RotData.RotMat[9*k+5] * MagData.my[l+N_cum] \
				   + RotData.RotMat[9*k+8] * MagData.mz[l+N_cum];

				// atomic position composition
				//X = MagData.RotMat[0] * xr \
                //  + MagData.RotMat[3] * yr \
                //  + MagData.RotMat[6] * zr;
				Y = MagData.RotMat[1] * xr \
				   + MagData.RotMat[4] * yr \
				   + MagData.RotMat[7] * zr;
				Z = MagData.RotMat[2] * xr \
				   + MagData.RotMat[5] * yr \
				   + MagData.RotMat[8] * zr;

				// phase function
				Psi = Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i];

				// cosine and sine values
				cos_val = cos(Psi);
				sin_val = sin(Psi);

				// cosine and sine summations
				nuc_real += NucData.Nuc[l+N_cum] * cos_val;
				nuc_imag -= NucData.Nuc[l+N_cum] * sin_val;

				mx_real += mxr * cos_val;
				mx_imag -= mxr * sin_val;
				my_real += myr * cos_val;
				my_imag -= myr * sin_val;
				mz_real += mzr * cos_val;
				mz_imag -= mzr * sin_val;

			}

			// rotations of the magnetization fourier components real part
			Mx_real = MagData.RotMat[0] * mx_real + MagData.RotMat[3] * my_real + MagData.RotMat[6] * mz_real;
			My_real = MagData.RotMat[1] * mx_real + MagData.RotMat[4] * my_real + MagData.RotMat[7] * mz_real;
			Mz_real = MagData.RotMat[2] * mx_real + MagData.RotMat[5] * my_real + MagData.RotMat[8] * mz_real;

			// rotations of the magnetization fourier components imaginary part
			Mx_imag = MagData.RotMat[0] * mx_imag + MagData.RotMat[3] * my_imag + MagData.RotMat[6] * mz_imag;
			My_imag = MagData.RotMat[1] * mx_imag + MagData.RotMat[4] * my_imag + MagData.RotMat[7] * mz_imag;
			Mz_imag = MagData.RotMat[2] * mx_imag + MagData.RotMat[5] * my_imag + MagData.RotMat[8] * mz_imag;


			cos_theta = cosf(SANSData.theta_2D[i]);
			sin_theta = sinf(SANSData.theta_2D[i]);

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
