#pragma once

// inline void FreeData(InputFileData* InputData,
//                      NuMagSANSData* Data){

// 	LogSystem::write("free memory...");

// 	if(InputData->NucData_activate_flag){
// 		free_NuclearData(&Data->NucData, &Data->NucData_gpu);
// 	}

// 	if(InputData->MagData_activate_flag){
// 		free_MagnetizationData(&Data->MagData, &Data->MagData_gpu);
// 	}

// 	if(InputData->StructData_activate_flag){
// 		free_StructureData(&Data->StructData, &Data->StructData_gpu);
// 	}

// 	if(InputData->RotData_activate_flag){
// 		free_RotationData(&Data->RotData, &Data->RotData_gpu);
// 	}

// 	free_ScatteringData(&Data->SANSData, &Data->SANSData_gpu);
// 	free_SpectralData(&Data->SpecData, &Data->SpecData_gpu);
// }


inline void FreeData(InputFileData* InputData,
					 MagnetizationData* MagData, MagnetizationData* MagData_gpu,
				     NuclearData* NucData, NuclearData* NucData_gpu,
				     StructureData* StructData, StructureData* StructData_gpu,
					 RotationData* RotData, RotationData* RotData_gpu,
					 ScatteringData* SANSData, ScatteringData* SANSData_gpu,
					 SpectralData* SpecData, SpectralData* SpecData_gpu){


	if(InputData->NucData_activate_flag){
		free_NuclearData(NucData, NucData_gpu);
	}

	if(InputData->MagData_activate_flag){
		free_MagnetizationData(MagData, MagData_gpu);
	}

	if(InputData->StructData_activate_flag){
		free_StructureData(StructData, StructData_gpu);
	}

	if(InputData->RotData_activate_flag){
		free_RotationData(RotData, RotData_gpu);
	}
	
	// free scattering data
	free_ScatteringData(SANSData, SANSData_gpu);

	// free spectral data
	free_SpectralData(SpecData, SpecData_gpu);



}