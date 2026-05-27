#pragma once

inline void FreeData(InputFileData* InputData,
                     NuMagSANSData* Data){

	LogSystem::write("free memory...");

	if(InputData->NucData_activate_flag){
		free_NuclearData(&Data->NucData, &Data->NucData_gpu);
	}

	if(InputData->MagData_activate_flag){
		free_MagnetizationData(&Data->MagData, &Data->MagData_gpu);
	}

	if(InputData->StructData_activate_flag){
		free_StructureData(&Data->StructData, &Data->StructData_gpu);
	}

	if(InputData->RotData_activate_flag){
		free_RotationData(&Data->RotData, &Data->RotData_gpu);
	}

	free_ScatteringData(&Data->SANSData, &Data->SANSData_gpu);
	free_SpectralData(&Data->SpecData, &Data->SpecData_gpu);
}
