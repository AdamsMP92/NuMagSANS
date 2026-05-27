#pragma once

inline void CheckKernelLaunch(cudaError_t err){
	if (err != cudaSuccess) {
		LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
	}
}

inline void SelectKernelRun(InputFileData* InputData,
                            NuMagSANSData* Data){

	cudaError_t err;
	NuclearData* NucData_gpu = &Data->NucData_gpu;
	MagnetizationData* MagData_gpu = &Data->MagData_gpu;
	StructureData* StructData_gpu = &Data->StructData_gpu;
	RotationData* RotData_gpu = &Data->RotData_gpu;
	ScatteringData* SANSData_gpu = &Data->SANSData_gpu;
	int L = (*SANSData_gpu->N_q) * (*SANSData_gpu->N_theta);
	LogSystem::write("total number of Fourier space bins: " + std::to_string(L));

	// Pure Magnetic Scattering Calculator without structure data #############################
	if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 0 && InputData->StructData_activate_flag == 0 && InputData->RotData_activate_flag == 0){
		LogSystem::write("run: Atomistic_MagSANS_Kernel_dilute");
		Atomistic_MagSANS_Kernel_dilute<<<(L+255)/256, 256>>>(*MagData_gpu, *SANSData_gpu);
		cudaDeviceSynchronize();
		err = cudaGetLastError();
		CheckKernelLaunch(err);
	}

	// Pure Nuclear Scattering Calculator without structure data ##############################
	if(InputData->MagData_activate_flag == 0 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 0 && InputData->RotData_activate_flag == 0){
		LogSystem::write("run: Atomistic_NucSANS_Kernel_dilute");
		Atomistic_NucSANS_Kernel_dilute<<<(L+255)/256, 256>>>(*NucData_gpu, *SANSData_gpu);
		cudaDeviceSynchronize();
		err = cudaGetLastError();
		CheckKernelLaunch(err);
	}

	// Combined Magnetic and Nuclear Scattering Calculator without structure data #############
	if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 0 && InputData->RotData_activate_flag == 0){
		LogSystem::write("run: Atomistic_NuMagSANS_Kernel_dilute");
		Atomistic_NuMagSANS_Kernel_dilute<<<(L+255)/256, 256>>>(*NucData_gpu, *MagData_gpu, *SANSData_gpu);
		cudaDeviceSynchronize();
		err = cudaGetLastError();
		CheckKernelLaunch(err);
	}

	// Pure Magnetic Scattering Calculator with structure data ################################
	if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 0 && InputData->StructData_activate_flag == 1 && InputData->RotData_activate_flag == 0){
		LogSystem::write("run: Atomistic_MagSANS_Kernel");
		Atomistic_MagSANS_Kernel<<<(L+255)/256, 256>>>(*MagData_gpu, *StructData_gpu, *SANSData_gpu);
		cudaDeviceSynchronize();
		err = cudaGetLastError();
		CheckKernelLaunch(err);
	}

	// Pure Nuclear Scattering Calculator with structure data #################################
	if(InputData->MagData_activate_flag == 0 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 1 && InputData->RotData_activate_flag == 0){
		LogSystem::write("run: Atomistic_NucSANS_Kernel");
		Atomistic_NucSANS_Kernel<<<(L+255)/256, 256>>>(*NucData_gpu, *StructData_gpu, *SANSData_gpu);
		cudaDeviceSynchronize();
		err = cudaGetLastError();
		CheckKernelLaunch(err);
	}

	// Combined Magnetic and Nuclear Scattering Calculator with structure data ################
	if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 1 && InputData->RotData_activate_flag == 0){
		LogSystem::write("run: Atomistic_NuMagSANS_Kernel");
		Atomistic_NuMagSANS_Kernel<<<(L+255)/256, 256>>>(*NucData_gpu, *MagData_gpu, *StructData_gpu, *SANSData_gpu);
		cudaDeviceSynchronize();
		err = cudaGetLastError();
		CheckKernelLaunch(err);
	}

	// Pure Magnetic Scattering Calculator with rotation data and without structure data ######
	if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 0 && InputData->StructData_activate_flag == 0 && InputData->RotData_activate_flag == 1){
		LogSystem::write("run: Atomistic_MagSANS_Kernel_RotDilute");
		Atomistic_MagSANS_Kernel_RotDilute<<<(L+255)/256, 256>>>(*MagData_gpu, *RotData_gpu, *SANSData_gpu);
		cudaDeviceSynchronize();
		err = cudaGetLastError();
		CheckKernelLaunch(err);
	}

	// Pure Nuclear Scattering Calculator with rotation data and without structure data #######
	if(InputData->MagData_activate_flag == 0 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 0 && InputData->RotData_activate_flag == 1){
		LogSystem::write("run: Atomistic_NucSANS_Kernel_RotDilute");
		Atomistic_NucSANS_Kernel_RotDilute<<<(L+255)/256, 256>>>(*NucData_gpu, *RotData_gpu, *SANSData_gpu);
		cudaDeviceSynchronize();
		err = cudaGetLastError();
		CheckKernelLaunch(err);
	}

	// Combined Magnetic and Nuclear Scattering Calculator with rotation data #################
	if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 0 && InputData->RotData_activate_flag == 1){
		LogSystem::write("run: Atomistic_NuMagSANS_Kernel_RotDilute");
		Atomistic_NuMagSANS_Kernel_RotDilute<<<(L+255)/256, 256>>>(*NucData_gpu, *MagData_gpu, *RotData_gpu, *SANSData_gpu);
		cudaDeviceSynchronize();
		err = cudaGetLastError();
		CheckKernelLaunch(err);
	}

	// Pure Magnetic Scattering Calculator with structure and rotation data ###################
	if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 0 && InputData->StructData_activate_flag == 1 && InputData->RotData_activate_flag == 1){
		LogSystem::write("run: Atomistic_MagSANS_Kernel_StructRot");
		Atomistic_MagSANS_Kernel_StructRot<<<(L+255)/256, 256>>>(*MagData_gpu, *StructData_gpu, *RotData_gpu, *SANSData_gpu);
		cudaDeviceSynchronize();
		err = cudaGetLastError();
		CheckKernelLaunch(err);
	}

	// Pure Nuclear Scattering Calculator with structure and rotation data ####################
	if(InputData->MagData_activate_flag == 0 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 1 && InputData->RotData_activate_flag == 1){
		LogSystem::write("run: Atomistic_NucSANS_Kernel_StructRot");
		Atomistic_NucSANS_Kernel_StructRot<<<(L+255)/256, 256>>>(*NucData_gpu, *StructData_gpu, *RotData_gpu, *SANSData_gpu);
		cudaDeviceSynchronize();
		err = cudaGetLastError();
		CheckKernelLaunch(err);
	}

	// Combined Magnetic and Nuclear Scattering Calculator with structure and rotation data ###
	if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 1 && InputData->RotData_activate_flag == 1){
		LogSystem::write("run: Atomistic_NuMagSANS_Kernel_StructRot");
		Atomistic_NuMagSANS_Kernel_StructRot<<<(L+255)/256, 256>>>(*NucData_gpu, *MagData_gpu, *StructData_gpu, *RotData_gpu, *SANSData_gpu);
		cudaDeviceSynchronize();
		err = cudaGetLastError();
		CheckKernelLaunch(err);
	}
}
