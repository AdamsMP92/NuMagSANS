#pragma once

// inline void InitializeData(InputFileData* InputData,
// 							   NucDataProperties* NucDataProp,
// 							   MagDataProperties* MagDataProp,
// 							   StructDataProperties* StructDataProp,
// 							   RotDataProperties* RotDataProp,
// 							   NuMagSANSData* Data,
// 							   int Data_File_Index){

// 		cudaError_t err;

// 		// initialize nuclear data ################################################################
// 		if(InputData->NucData_activate_flag){
// 			init_NuclearData(&Data->NucData, &Data->NucData_gpu, NucDataProp, InputData, Data_File_Index);
// 			cudaDeviceSynchronize();
// 			err = cudaGetLastError();
// 			CheckKernelLaunch(err);
// 		}
		
// 		// initialize magnetization data ##########################################################
// 		if(InputData->MagData_activate_flag){
// 			init_MagnetizationData(&Data->MagData, &Data->MagData_gpu, MagDataProp, InputData, Data_File_Index);
// 			cudaDeviceSynchronize();
// 			err = cudaGetLastError();
// 			CheckKernelLaunch(err);
// 			//disp_MagnetizationData(&Data->MagData);
// 		}

// 		// initialize structure data ##############################################################
// 		if(InputData->StructData_activate_flag){
// 			init_StructureData(&Data->StructData, &Data->StructData_gpu, StructDataProp, InputData);
// 			cudaDeviceSynchronize();
// 			err = cudaGetLastError();
// 			CheckKernelLaunch(err);
// 			//disp_StructureData(&Data->StructData);
// 		}

// 		// initialize rotation data ###############################################################
// 		if(InputData->RotData_activate_flag){
// 			init_RotationData(&Data->RotData, &Data->RotData_gpu, RotDataProp, InputData);
// 			cudaDeviceSynchronize();
// 			err = cudaGetLastError();
// 			CheckKernelLaunch(err);
// 		}

// 		// initialize scattering data #############################################################
// 		init_ScatteringData(InputData, &Data->SANSData, &Data->SANSData_gpu);
// 		cudaDeviceSynchronize();
// 		err = cudaGetLastError();
// 		CheckKernelLaunch(err);

// 		// initialize spectral data ###############################################################
// 		init_SpectralData(InputData, &Data->SANSData, &Data->SpecData, &Data->SpecData_gpu);
// 		cudaDeviceSynchronize();
// 		err = cudaGetLastError();
// 		CheckKernelLaunch(err);
// 	}

inline void InitializeData(InputFileData* InputData,
						   NucDataProperties* NucDataProp,
						   MagDataProperties* MagDataProp,
						   StructDataProperties* StructDataProp,
						   RotDataProperties* RotDataProp,
						   NuclearData* NucData, NuclearData* NucData_gpu,
						   MagnetizationData* MagData, MagnetizationData* MagData_gpu,
						   StructureData* StructData, StructureData* StructData_gpu,
						   RotationData* RotData, RotationData* RotData_gpu,
						   ScatteringData* SANSData, ScatteringData* SANSData_gpu,
						   SpectralData* SpecData, SpectralData* SpecData_gpu,
						   int Data_File_Index){

	cudaError_t err;

	// initialize nuclear data ################################################################
	if(InputData->NucData_activate_flag){
		init_NuclearData(NucData, NucData_gpu, NucDataProp, InputData, Data_File_Index);
		cudaDeviceSynchronize();
		err = cudaGetLastError();
		if (err != cudaSuccess) {
			LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
		}
	}
	
	// initialize magnetization data ##########################################################
	if(InputData->MagData_activate_flag){
		init_MagnetizationData(MagData, MagData_gpu, MagDataProp, InputData, Data_File_Index);
		cudaDeviceSynchronize();
		err = cudaGetLastError();
		if (err != cudaSuccess) {
			LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
		}
		//disp_MagnetizationData(MagData);
	}

	// initialize structure data ##############################################################
	if(InputData->StructData_activate_flag){
		init_StructureData(StructData, StructData_gpu, StructDataProp, InputData);
		cudaDeviceSynchronize();
		err = cudaGetLastError();
		if (err != cudaSuccess) {
			LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
		}
		//disp_StructureData(StructData);
	}

	// initialize rotation data ###############################################################
	if(InputData->RotData_activate_flag){
		init_RotationData(RotData, RotData_gpu, RotDataProp, InputData);
		cudaDeviceSynchronize();
		err = cudaGetLastError();
		if (err != cudaSuccess) {
			LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
		}
	}

	// initialize scattering data #############################################################
	init_ScatteringData(InputData, SANSData, SANSData_gpu);
	cudaDeviceSynchronize();
	err = cudaGetLastError();
	if (err != cudaSuccess) {
		LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
	}

	// initialize spectral data ###############################################################
	init_SpectralData(InputData, SANSData, SpecData, SpecData_gpu);
	cudaDeviceSynchronize();
	err = cudaGetLastError();
	if (err != cudaSuccess) {
		LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
	}
}
