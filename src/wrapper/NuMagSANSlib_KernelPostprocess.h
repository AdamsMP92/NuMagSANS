#pragma once

inline void KernelPostprocessRunMulti(InputFileData* InputData,
									  ScatteringData* SANSData,
									  ScatteringData* SANSData_gpu,
									  SpectralData* SpecData_gpu){

	int L = (*SANSData->N_q) * (*SANSData->N_theta);

	bool compute_1D_azimuthal_average = any_active(InputData->OutFlags.SANS1D);
	bool compute_1D_corr = any_active(InputData->OutFlags.Corr1D);
	bool compute_1D_pair = any_active(InputData->OutFlags.PairDist1D);

	if(compute_1D_azimuthal_average || compute_1D_corr || compute_1D_pair){
		LogSystem::write("run: azimuthal averaging");
		AzimuthalAverage<<<(L+255)/256, 256>>>(*SANSData_gpu);
		CheckCUDAKernelRun("AzimuthalAverage");
	}

	if(compute_1D_corr || compute_1D_pair){
		LogSystem::write("run: 1D correlation functions");
		DistributionFunctions<<<(L+255), 256>>>(*SANSData_gpu);
		CheckCUDAKernelRun("DistributionFunctions");
	}

	bool compute_2D_correlation = any_active(InputData->OutFlags.Corr2D);

	if(compute_2D_correlation){
		LogSystem::write("run: 2D correlation functions");
		CorrelationFunction_2D<<<(L+255)/256, 256>>>(*SANSData_gpu);
		CheckCUDAKernelRun("CorrelationFunction_2D");
	}

	if(InputData->AngularSpec_activate_flag){
		LogSystem::write("run: angular spectrum analyzer");
		ComputeSpectralDecomposition<<<(L+255)/256, 256>>>(*SANSData_gpu, *SpecData_gpu);
		CheckCUDAKernelRun("ComputeSpectralDecomposition");

		LogSystem::write("run: angular amplitude spectrum analyzer");
		ComputeAngularSpectrumAmplitudes<<<(L+255)/256, 256>>>(*SpecData_gpu);
		CheckCUDAKernelRun("ComputeAngularSpectrumAmplitudes");
	}
}


inline void KernelPostprocessRun(InputFileData* InputData,
								 NuMagSANSData* Data){

	KernelPostprocessRunMulti(InputData,
							  &Data->SANSData,
							  &Data->SANSData_gpu,
							  &Data->SpecData_gpu);

}
