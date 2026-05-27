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
#include "NuMagSANSlib_MemoryInfo.h"
#include "NuMagSANSlib_TimeMeasure.h"
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
#include "NuMagSANSlib_gpuKernel.h"
#include "NuMagSANSlib_InitializeData.h"
#include "NuMagSANSlib_KernelSelector.h"
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

	cudaError_t err;
	
	LogSystem::write("################################################################################");
	LogSystem::write("## Run - NuMagSANS #############################################################");
	LogSystem::write("################################################################################");
	LogSystem::write("");

	GPUMemoryInfo MemoryBeforeRun = GetGPUMemoryInfo();
 	LogCurrentGPUMemoryDifference(MemoryBeforeRun);

	// start time measurement #################################################################
	TimeMeasure TotalTime = StartTimeMeasure();


	NuclearData NucData, NucData_gpu;
	MagnetizationData MagData, MagData_gpu;
    StructureData StructData, StructData_gpu;
	RotationData RotData, RotData_gpu;
	ScatteringData SANSData, SANSData_gpu;
	SpectralData SpecData, SpecData_gpu;

	InitializeData(InputData,
				   NucDataProp,
				   MagDataProp,
				   StructDataProp,
				   RotDataProp,
				   &NucData, &NucData_gpu,
				   &MagData, &MagData_gpu,
				   &StructData, &StructData_gpu,
				   &RotData, &RotData_gpu,
				   &SANSData, &SANSData_gpu,
				   &SpecData, &SpecData_gpu,
				   Data_File_Index);

	LogCurrentGPUMemoryDifference(MemoryBeforeRun);
	
	// initialize scaling factors #############################################################
	ScalingFactors ScalFactors;
	init_ScalingFactors(&ScalFactors, InputData, &MagData, &NucData, &SANSData);
	
	SelectKernelRun(InputData,
					&NucData_gpu,
					&MagData_gpu,
					&StructData_gpu,
					&RotData_gpu,
					&SANSData,
					&SANSData_gpu);

	// // compute azimuthal average 1D ###########################################################

	// int L = (*SANSData.N_q) * (*SANSData.N_theta);
	// bool compute_1D_azimuthal_average = any_active(InputData->OutFlags.SANS1D);
	// bool compute_1D_corr = any_active(InputData->OutFlags.Corr1D);
	// bool compute_1D_pair = any_active(InputData->OutFlags.PairDist1D);

	// if(compute_1D_azimuthal_average || compute_1D_corr || compute_1D_pair){

	// 	LogSystem::write("run: azimuthal averaging");
	// 	AzimuthalAverage<<<(L+255)/256, 256>>>(SANSData_gpu);
	// 	cudaDeviceSynchronize();
	// 	err = cudaGetLastError();
	// 	if (err != cudaSuccess) {
	// 	   LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err)); 
	// 	}

	// }

	// // compute 1D pair distance distribution and correlation function #########################
	// if(compute_1D_corr || compute_1D_pair){
	// 	LogSystem::write("run: 1D correlation functions");
	// 	DistributionFunctions<<<(L+255), 256>>>(SANSData_gpu);
	// 	cudaDeviceSynchronize();
	// 	err = cudaGetLastError();
	// 	if (err != cudaSuccess) {
	//    		LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
	// 	}
	// }

	// // compute 2D correlation functions ######################################################
	// bool compute_2D_correlation = any_active(InputData->OutFlags.Corr2D);	

	// if(compute_2D_correlation){
	// 	LogSystem::write("run: 2D correlation functions");
	// 	CorrelationFunction_2D<<<(L+255)/256, 256>>>(SANSData_gpu);
	// 	cudaDeviceSynchronize();
	// 	err = cudaGetLastError();
	// 	if (err != cudaSuccess) {
	// 	   LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
	// 	}
	// }


	// // compute angular spectral intensities ###################################################
	// if(InputData->AngularSpec_activate_flag){
	
	// 	LogSystem::write("run: angular spectrum analyzer");
	// 	ComputeSpectralDecomposition<<<(L+255)/256, 256>>>(SANSData_gpu, SpecData_gpu);
	// 	err = cudaGetLastError();
	// 	if (err != cudaSuccess) {
	// 	   LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
	// 	}

	// // compute angular spectral amplitudes ####################################################
    //    	LogSystem::write("run: angular amplitude spectrum analyzer");
	// 	ComputeAngularSpectrumAmplitudes<<<(L+255)/256, 256>>>(SpecData_gpu);
	// 	err = cudaGetLastError();
    //    	if (err != cudaSuccess) {
    //       		LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
    //    	}

	// }

	KernelPostprocessRun(InputData,
						 &SANSData,
						 &SANSData_gpu,
						 &SpecData_gpu);


	// copy scattering data from GPU to RAM ###################################################
	ExportData(InputData,
			   &ScalFactors,
			   &SANSData, &SANSData_gpu,
			   &SpecData, &SpecData_gpu,
			   Data_File_Index);

	// free memory ############################################################################
	FreeData(InputData,
			 &MagData, &MagData_gpu,
			 &NucData, &NucData_gpu,
			 &StructData, &StructData_gpu,
			 &RotData, &RotData_gpu,
			 &SANSData, &SANSData_gpu,
			 &SpecData, &SpecData_gpu);
	
	// print result of time measurement #######################################################
    LogElapsedTime(TotalTime);

	LogCurrentGPUMemoryDifference(MemoryBeforeRun);

	LogSystem::write("################################################################################");
	LogSystem::write("## Finished - NuMagSANS ########################################################");
	LogSystem::write("################################################################################");
	LogSystem::write("");
	LogSystem::write("");
}






// // File         : NuMagSANSlib.h
// // Author       : Dr. Michael Philipp ADAMS 
// // Company      : University of Luxembourg
// // Department   : Department of Physics and Materials Sciences
// // Group        : NanoMagnetism Group
// // Group Leader : Prof. Andreas Michels
// // Version      : 29 October 2025
// // OS           : Linux Ubuntu
// // Language     : CUDA C++

// #include <iostream>
// #include <fstream>
// #include <sstream>
// #include <sys/stat.h>
// #include <sys/types.h>
// #include <math.h>
// #include <string>
// #include <vector>
// #include <stdlib.h>
// #include <time.h>
// #include <cuda_runtime.h>
// #include <cublas_v2.h>
// #include <stdexcept>
// #include <math.h>
// #include <chrono>
// #include <dirent.h>
// #include <unistd.h>

// //#pragma once
// #include "NuMagSANSlib_LogFile.h"
// #include "NuMagSANSlib_HelperFun.h"
// #include "NuMagSANSlib_StringCompare.h"
// #include "NuMagSANSlib_ReadWrite.h"
// #include "NuMagSANSlib_Directory.h"
// #include "NuMagSANSlib_InputFileInterpreter.h"
// #include "NuMagSANSlib_MagDataExplorer.h"
// #include "NuMagSANSlib_NucDataExplorer.h"
// #include "NuMagSANSlib_StructureDataExplorer.h"
// #include "NuMagSANSlib_RotationDataExplorer.h"
// #include "NuMagSANSlib_MagData.h"
// #include "NuMagSANSlib_NucData.h"
// #include "NuMagSANSlib_StructureData.h"
// #include "NuMagSANSlib_RotationData.h"
// #include "NuMagSANSlib_SANSData.h"
// #include "NuMagSANSlib_SpectralData.h"
// #include "NuMagSANSlib_Data.h"
// #include "NuMagSANSlib_MemoryInfo.h"
// #include "NuMagSANSlib_TimeMeasure.h"
// #include "NuMagSANSlib_gpuKernel.h"
// #include "NuMagSANSlib_KernelSelector.h"
// #include "NuMagSANSlib_InitializeData.h"
// #include "NuMagSANSlib_KernelPostprocess.h"
// #include "NuMagSANSlib_ExportData.h"
// #include "NuMagSANSlib_FreeData.h"

// using namespace std;


// void NuMagSANS_Calculator(InputFileData* InputData, \
// 					      NucDataProperties* NucDataProp,\
//                           MagDataProperties* MagDataProp,\
//                           StructDataProperties* StructDataProp, \
// 						  RotDataProperties* RotDataProp, \
//                           int Data_File_Index){
	
// 	LogSystem::write("################################################################################");
// 	LogSystem::write("## Run - NuMagSANS #############################################################");
// 	LogSystem::write("################################################################################");
// 	LogSystem::write("");

// 	// initial gpu memory status ######################################
// 	GPUMemoryInfo MemoryBeforeRun = GetGPUMemoryInfo();
// 	LogCurrentGPUMemoryDifference(MemoryBeforeRun);

// 	// start time measurement #################################################################
// 	TimeMeasure TotalTime = StartTimeMeasure();

// 	// define host and GPU data objects #######################################################
// 	NuMagSANSData Data = {};

// 	// initialize host and GPU data ###########################################################
// 	InitializeData(InputData, NucDataProp, MagDataProp, StructDataProp, RotDataProp,
// 				   &Data, Data_File_Index);

// 	// report GPU memory state after data initialization #######################################
// 	LogCurrentGPUMemoryDifference(MemoryBeforeRun);
	
// 	// initialize scaling factors #############################################################
// 	ScalingFactors ScalFactors;
// 	init_ScalingFactors(&ScalFactors, InputData, &Data.MagData, &Data.NucData, &Data.SANSData);
	
// 	// compute 2D SANS cross sections #########################################################
// 	SelectKernelRun(InputData,
// 					Data.NucData_gpu,
// 					Data.MagData_gpu,
// 					Data.StructData_gpu,
// 					Data.RotData_gpu,
// 					Data.SANSData_gpu);
	
// 	// run postprocessing kernels #############################################################
// 	KernelPostprocessRun(InputData, Data.SANSData_gpu, Data.SpecData_gpu);

// 	// export scattering and spectral data ####################################################
// 	ExportData(InputData, &ScalFactors, &Data, Data_File_Index);

// 	// free memory ############################################################################
// 	FreeData(InputData, &Data);
	
// 	// print result of time measurement #######################################################
// 	LogElapsedTime(TotalTime);

// 	LogCurrentGPUMemoryDifference(MemoryBeforeRun);

// 	LogSystem::write("################################################################################");
// 	LogSystem::write("## Finished - NuMagSANS ########################################################");
// 	LogSystem::write("################################################################################");
// 	LogSystem::write("");

// }
