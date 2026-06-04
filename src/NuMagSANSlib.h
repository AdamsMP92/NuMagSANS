// File         : NuMagSANSlib.h
// Author       : Dr. Michael Philipp ADAMS
// Company      : University of Luxembourg
// Department   : Department of Physics and Materials Sciences
// Group        : NanoMagnetism Group
// Group Leader : Prof. Andreas Michels
// Version      : 29 October 2025
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

//#pragma once
#include "helper/NuMagSANSlib_LogFile.h"
#include "helper/NuMagSANSlib_MemoryInfo.h"
#include "helper/NuMagSANSlib_TimeMeasure.h"
#include "helper/NuMagSANSlib_CUDAError.h"
#include "helper/NuMagSANSlib_HelperFun.h"
#include "helper/NuMagSANSlib_StringCompare.h"
#include "helper/NuMagSANSlib_ReadWrite.h"
#include "helper/NuMagSANSlib_Directory.h"
#include "InputData/NuMagSANSlib_InputFileInterpreter.h"
#include "InputData/MagData/NuMagSANSlib_MagDataExplorer.h"
#include "InputData/NucData/NuMagSANSlib_NucDataExplorer.h"
#include "InputData/StructData/NuMagSANSlib_StructureDataExplorer.h"
#include "InputData/RotData/NuMagSANSlib_RotationDataExplorer.h"
#include "InputData/MagData/NuMagSANSlib_MagData.h"
#include "InputData/NucData/NuMagSANSlib_NucData.h"
#include "InputData/StructData/NuMagSANSlib_StructureData.h"
#include "InputData/RotData/NuMagSANSlib_RotationData.h"
#include "OutputData/SANSData/NuMagSANSlib_SANSData.h"
#include "OutputData/SpectralData/NuMagSANSlib_SpectralData.h"
#include "wrapper/NuMagSANSlib_Data.h"
#include "gpu_kernels/NuMagSANSlib_gpuKernel.h"
#include "wrapper/NuMagSANSlib_InitializeData.h"
#include "wrapper/NuMagSANSlib_KernelSelector.h"
#include "wrapper/NuMagSANSlib_KernelPostprocess.h"
#include "wrapper/NuMagSANSlib_ExportData.h"
#include "wrapper/NuMagSANSlib_FreeData.h"


using namespace std;

void NuMagSANS_Calculator(InputFileData* InputData, \
					      NucDataProperties* NucDataProp,\
                          MagDataProperties* MagDataProp,\
                          StructDataProperties* StructDataProp, \
						  RotDataProperties* RotDataProp, \
                          int Data_File_Index){

	LogSystem::write("################################################################################");
	LogSystem::write("## Run - NuMagSANS #############################################################");
	LogSystem::write("################################################################################");
	LogSystem::write("");

	GPUMemoryInfo MemoryBeforeRun = GetGPUMemoryInfo();
 	LogCurrentGPUMemoryDifference(MemoryBeforeRun);

	// start time measurement
	TimeMeasure TotalTime = StartTimeMeasure();

	// initialize data
	NuMagSANSData Data;
	InitializeData(InputData,
		           NucDataProp, MagDataProp, StructDataProp, RotDataProp,
				   &Data, Data_File_Index);

	// check memory
	LogCurrentGPUMemoryDifference(MemoryBeforeRun);

	if(InputData->RotData_activate_flag && InputData->RotDataLoop_flag){

		for(int RotDataLoop_Counter = 0;
			RotDataLoop_Counter < InputData->RotDataLoop_IndexArray.size();
			RotDataLoop_Counter++){

			int RotData_File_Index = InputData->RotDataLoop_IndexArray[RotDataLoop_Counter];

			LogSystem::write("inner RotData loop active: RotData index "
			                 + std::to_string(RotData_File_Index));

			if(RotDataLoop_Counter != 0){
				ResetSANSOutputData(InputData, &Data);
			}

			SetActiveRotDataFile(RotDataProp, RotData_File_Index);
			new_read_RotationData(&Data.RotData, &Data.RotData_gpu, RotDataProp, InputData);
			CheckCUDALastError("reload rotation data");

			// run gpu kernel
			SelectKernelRun(InputData, &Data);

			// run gpu post-processing kernels
			KernelPostprocessRun(InputData, &Data);

			// export data
			ExportData(InputData, &Data, Data_File_Index, RotData_File_Index);
			CheckCUDALastError("export scattering data");
		}

	}else{

		// run gpu kernel
		SelectKernelRun(InputData, &Data);

		// run gpu post-processing kernels
		KernelPostprocessRun(InputData, &Data);

		// export data
		ExportData(InputData, &Data, Data_File_Index);
		CheckCUDALastError("export scattering data");
	}

	// free memory
	FreeData(InputData, &Data);
	CheckCUDALastError("free data");

	// print result of time measurement
    LogElapsedTime(TotalTime);

	LogCurrentGPUMemoryDifference(MemoryBeforeRun);

	LogSystem::write("################################################################################");
	LogSystem::write("## Finished - NuMagSANS ########################################################");
	LogSystem::write("################################################################################");
	LogSystem::write("");
	LogSystem::write("");
}
