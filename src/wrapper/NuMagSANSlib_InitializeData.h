#pragma once

inline void InitializeDataMulti(InputFileData* InputData,
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
									ScalingFactors* ScalFactors,
									int Data_File_Index){

	// initialize nuclear data ################################################################
	if(InputData->NucData_activate_flag){
		init_NuclearData(NucData, NucData_gpu, NucDataProp, InputData, Data_File_Index);
		CheckCUDALastError("initialize nuclear data");
	}
	
	// initialize magnetization data ##########################################################
	if(InputData->MagData_activate_flag){
		init_MagnetizationData(MagData, MagData_gpu, MagDataProp, InputData, Data_File_Index);
		CheckCUDALastError("initialize magnetization data");
		//disp_MagnetizationData(MagData);
	}

	// initialize structure data ##############################################################
	if(InputData->StructData_activate_flag){
		if(InputData->StructDataLoop_flag){
			init_StructureDataMemory(StructData, StructData_gpu, StructDataProp, InputData);
		}else{
			init_StructureData(StructData, StructData_gpu, StructDataProp, InputData);
		}
		CheckCUDALastError("initialize structure data");
		//disp_StructureData(StructData);
	}

	// initialize rotation data ###############################################################
	if(InputData->RotData_activate_flag){
		if(InputData->RotDataLoop_flag){
			init_RotationDataMemory(RotData, RotData_gpu, RotDataProp, InputData);
		}else{
			init_RotationData(RotData, RotData_gpu, RotDataProp, InputData);
		}
		CheckCUDALastError("initialize rotation data");
	}

	// initialize scattering data #############################################################
	init_ScatteringData(InputData, SANSData, SANSData_gpu);
	CheckCUDALastError("initialize scattering data");

	// initialize spectral data ###############################################################
	init_SpectralData(InputData, SANSData, SpecData, SpecData_gpu);
	CheckCUDALastError("initialize spectral data");

	// initialize scaling factors
	init_ScalingFactors(ScalFactors, InputData, MagData, NucData, SANSData);

}



inline void InitializeData(InputFileData* InputData,
							   NucDataProperties* NucDataProp,
							   MagDataProperties* MagDataProp,
							   StructDataProperties* StructDataProp,
							   RotDataProperties* RotDataProp, 
							   NuMagSANSData* Data,
							   int Data_File_Index){

		InitializeDataMulti(InputData,
							NucDataProp,
							MagDataProp,
							StructDataProp,
							RotDataProp,
							&Data->NucData, &Data->NucData_gpu,
							&Data->MagData, &Data->MagData_gpu,
							&Data->StructData, &Data->StructData_gpu,
							&Data->RotData, &Data->RotData_gpu,
							&Data->SANSData, &Data->SANSData_gpu,
							&Data->SpecData, &Data->SpecData_gpu,
							&Data->ScalFactors,
							Data_File_Index);


	}

inline void ResetSANSOutputData(InputFileData* InputData,
								NuMagSANSData* Data){

	free_ScatteringData(&Data->SANSData, &Data->SANSData_gpu);
	free_SpectralData(&Data->SpecData, &Data->SpecData_gpu);

	init_ScatteringData(InputData, &Data->SANSData, &Data->SANSData_gpu);
	init_SpectralData(InputData, &Data->SANSData, &Data->SpecData, &Data->SpecData_gpu);

	init_ScalingFactors(&Data->ScalFactors, InputData, &Data->MagData, &Data->NucData, &Data->SANSData);

	CheckCUDALastError("reset SANS output data");
}
