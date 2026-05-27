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
#include "NuMagSANSlib_LogFile.h"
#include "NuMagSANSlib_HelperFun.h"
#include "NuMagSANSlib_StringCompare.h"
#include "NuMagSANSlib_ReadWrite.h"
#include "NuMagSANSlib_Directory.h"
#include "NuMagSANSlib_InputFileInterpreter.h"
#include "NuMagSANSlib_MagDataExplorer.h"
#include "NuMagSANSlib_NucDataExplorer.h"
#include "NuMagSANSlib_StructureDataExplorer.h"
#include "NuMagSANSlib_RotationDataExplorer.h"
#include "NuMagSANSlib_MagData.h"
#include "NuMagSANSlib_NucData.h"
#include "NuMagSANSlib_StructureData.h"
#include "NuMagSANSlib_RotationData.h"
#include "NuMagSANSlib_SANSData.h"
#include "NuMagSANSlib_SpectralData.h"
#include "NuMagSANSlib_Data.h"
#include "NuMagSANSlib_MemoryInfo.h"
#include "NuMagSANSlib_TimeMeasure.h"
#include "NuMagSANSlib_gpuKernel.h"
#include "NuMagSANSlib_KernelSelector.h"
#include "NuMagSANSlib_InitializeData.h"
#include "NuMagSANSlib_KernelPostprocess.h"
#include "NuMagSANSlib_ExportData.h"
#include "NuMagSANSlib_FreeData.h"

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

	// initial gpu memory status ######################################
	GPUMemoryInfo MemoryBeforeRun = GetGPUMemoryInfo();
	LogCurrentGPUMemoryDifference(MemoryBeforeRun);

	// start time measurement #################################################################
	TimeMeasure TotalTime = StartTimeMeasure();

	// define host and GPU data objects #######################################################
	NuMagSANSData Data = {};

	// initialize host and GPU data ###########################################################
	InitializeData(InputData, NucDataProp, MagDataProp, StructDataProp, RotDataProp,
				   &Data, Data_File_Index);

	// report GPU memory state after data initialization #######################################
	LogCurrentGPUMemoryDifference(MemoryBeforeRun);
	
	// initialize scaling factors #############################################################
	ScalingFactors ScalFactors;
	init_ScalingFactors(&ScalFactors, InputData, &Data.MagData, &Data.NucData, &Data.SANSData);
	
	// compute 2D SANS cross sections #########################################################
	SelectKernelRun(InputData, &Data);
	
	// run postprocessing kernels #############################################################
	KernelPostprocessRun(InputData, &Data);

	// export scattering and spectral data ####################################################
	ExportData(InputData, &ScalFactors, &Data, Data_File_Index);

	// free memory ############################################################################
	FreeData(InputData, &Data);
	
	// print result of time measurement #######################################################
	LogElapsedTime(TotalTime);

	LogCurrentGPUMemoryDifference(MemoryBeforeRun);

	LogSystem::write("################################################################################");
	LogSystem::write("## Finished - NuMagSANS ########################################################");
	LogSystem::write("################################################################################");
	LogSystem::write("");

}
