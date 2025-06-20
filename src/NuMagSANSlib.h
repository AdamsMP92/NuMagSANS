// File         : NuMagSANSlib.h
// Author       : Michael Philipp ADAMS, M.Sc. 
// Company      : University of Luxembourg
// Department   : Department of Physics and Materials Sciences
// Group        : NanoMagnetism Group
// Group Leader : Prof. Andreas Michels
// Version      : 23 November 2024
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
#include "NuMagSANSlib_HelperFun.h"
#include "NuMagSANSlib_StringCompare.h"
#include "NuMagSANSlib_ReadWrite.h"
#include "NuMagSANSlib_Directory.h"
#include "NuMagSANSlib_InputFileInterpreter.h"
#include "NuMagSANSlib_MagDataExplorer.h"
#include "NuMagSANSlib_NucDataExplorer.h"
#include "NuMagSANSlib_StructureDataExplorer.h"
#include "NuMagSANSlib_MagData.h"
#include "NuMagSANSlib_NucData.h"
#include "NuMagSANSlib_StructureData.h"
#include "NuMagSANSlib_SANSData.h"
#include "NuMagSANSlib_gpuKernel.h"

using namespace std;

void NuMagSANS_Calculator(InputFileData* InputData, \
					      NucDataProperties* NucDataProp,\
                          MagDataProperties* MagDataProp,\
                          StructDataProperties* StructDataProp, \
                          int Data_File_Index){

	cout << "################################################################################" << "\n";
	cout << "## Run - NuMagSANS #############################################################" << "\n";
	cout << "################################################################################" << "\n\n";

	// start time measurement #################################################################
	auto start_total_time = std::chrono::high_resolution_clock::now();

	// initialize nuclear data ################################################################
	NuclearData NucData, NucData_gpu;
	if(InputData->NucData_activate_flag){
		init_NuclearData(&NucData, &NucData_gpu, NucDataProp, InputData, Data_File_Index);
	}
	
	// initialize magnetization data ##########################################################
	MagnetizationData MagData, MagData_gpu;
	if(InputData->MagData_activate_flag){
		init_MagnetizationData(&MagData, &MagData_gpu, MagDataProp, InputData, Data_File_Index);
		//disp_MagnetizationData(&MagData);
	}

	// initialize structure data ##############################################################
	StructureData StructData, StructData_gpu;
	if(InputData->StructData_activate_flag){
		init_StructureData(&StructData, &StructData_gpu, StructDataProp, InputData);
		//disp_StructureData(&StructData);
	}

	// initialize scattering data #############################################################
	ScatteringData SANSData, SANSData_gpu;
	init_ScatteringData(InputData, &SANSData, &SANSData_gpu);

	// initialize scaling factors #############################################################
	ScalingFactors ScalFactors;
	init_ScalingFactors(&ScalFactors, InputData, &MagData, &NucData, &SANSData);

	// compute 2D SANS cross sections #########################################################
	int L = (*SANSData.N_q) * (*SANSData.N_theta);
	cout << "L = " << L << "\n";
	cudaError_t err;

		// Pure Magnetic Scattering Calculator without structure data #####################################################################
		if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 0 && InputData->StructData_activate_flag == 0){
			cout << "Run: Atomistic_MagSANS_Kernel_dilute" << "\n";
			Atomistic_MagSANS_Kernel_dilute<<<(L+255)/256, 256>>>(MagData_gpu, SANSData_gpu);
			cudaDeviceSynchronize();
			err = cudaGetLastError();
			if (err != cudaSuccess) {
			std::cout << "Kernel launch failed: " << cudaGetErrorString(err) << std::endl;
			}
		}

		// Pure Nuclear Scattering Calculator without structure data #######################################################################
		if(InputData->MagData_activate_flag == 0 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 0){
			cout << "Run: Atomistic_NucSANS_Kernel_dilute" << "\n";
			Atomistic_NucSANS_Kernel_dilute<<<(L+255)/256, 256>>>(NucData_gpu, SANSData_gpu);
			cudaDeviceSynchronize();
			err = cudaGetLastError();
			if (err != cudaSuccess) {
			std::cout << "Kernel launch failed: " << cudaGetErrorString(err) << std::endl;
			}
		}

		// Combined Magnetic and Nuclear Scattering Calculator without structure data ######################################################
		if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 0){
			cout << "Run: Atomistic_NuMagSANS_Kernel_dilute" << "\n";
			Atomistic_NuMagSANS_Kernel_dilute<<<(L+255)/256, 256>>>(NucData_gpu, MagData_gpu, SANSData_gpu);
			cudaDeviceSynchronize();
			err = cudaGetLastError();
			if (err != cudaSuccess) {
			std::cout << "Kernel launch failed: " << cudaGetErrorString(err) << std::endl;
			}
		}


		// Pure Magnetic Scattering Calculator with structure data #########################################################################
		if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 0 && InputData->StructData_activate_flag == 1){
			cout << "Run: Atomistic_MagSANS_Kernel" << "\n";
			Atomistic_MagSANS_Kernel<<<(L+255)/256, 256>>>(MagData_gpu, StructData_gpu, SANSData_gpu);
			cudaDeviceSynchronize();
			err = cudaGetLastError();
			if (err != cudaSuccess) {
			std::cout << "Kernel launch failed: " << cudaGetErrorString(err) << std::endl;
			}
		}

		// Pure Nuclear Scattering Calculator with structure data ###########################################################################
		if(InputData->MagData_activate_flag == 0 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 1){
			cout << "Run: Atomistic_NucSANS_Kernel" << "\n";
			Atomistic_NucSANS_Kernel<<<(L+255)/256, 256>>>(NucData_gpu, StructData_gpu, SANSData_gpu);
			cudaDeviceSynchronize();
			err = cudaGetLastError();
			if (err != cudaSuccess) {
			std::cout << "Kernel launch failed: " << cudaGetErrorString(err) << std::endl;
			}
		}

		// Combined Magnetic and Nuclear Scattering Calculator with structure data #########################################################
		if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 1){
			cout << "Run: Atomistic_NuMagSANS_Kernel" << "\n";
			Atomistic_NuMagSANS_Kernel<<<(L+255)/256, 256>>>(NucData_gpu, MagData_gpu, StructData_gpu, SANSData_gpu);
			cudaDeviceSynchronize();
			err = cudaGetLastError();
			if (err != cudaSuccess) {
			std::cout << "Kernel launch failed: " << cudaGetErrorString(err) << std::endl;
			}
		}



	// compute azimuthal average 1D ###########################################################
	AzimuthalAverage<<<(L+255)/256, 256>>>(SANSData_gpu);
	cudaDeviceSynchronize();
	err = cudaGetLastError();
	if (err != cudaSuccess) {
	    std::cout << "Kernel launch failed: " << cudaGetErrorString(err) << std::endl;
	}

	// compute 1D pair distance distribution and correlation function #########################
	DistributionFunctions<<<(L+255), 256>>>(SANSData_gpu);
	cudaDeviceSynchronize();
	err = cudaGetLastError();
	if (err != cudaSuccess) {
	    std::cout << "Kernel launch failed: " << cudaGetErrorString(err) << std::endl;
	}

	// compute 2D correlation functions ######################################################
	CorrelationFunction_2D<<<(L+255)/256, 256>>>(SANSData_gpu);
	cudaDeviceSynchronize();
	err = cudaGetLastError();
	if (err != cudaSuccess) {
	    std::cout << "Kernel launch failed: " << cudaGetErrorString(err) << std::endl;
	}

	// copy scattering data from GPU to RAM ###################################################
	copyGPU2RAM_ScatteringData(&SANSData, &SANSData_gpu);

	// scaling of the scattering data on RAM ##################################################
	cout << "scaling of the scattering data..." << "\n";
	scale_ScatteringData(&ScalFactors, &SANSData, InputData);

	// write scattering data to csv files #####################################################
	//write2CSV_ScatteringData(InputData, &SANSData, Data_File_Index);
	write2CSVtable_ScatteringData(InputData, &SANSData, Data_File_Index);


	// free memory ############################################################################
	if(InputData->NucData_activate_flag){
		free_NuclearData(&NucData, &NucData_gpu);
	}

	if(InputData->MagData_activate_flag){
		free_MagnetizationData(&MagData, &MagData_gpu);
	}

	if(InputData->StructData_activate_flag){
		free_StructureData(&StructData, &StructData_gpu);
	}
	
	free_ScatteringData(&SANSData, &SANSData_gpu);

	// print result of time measurement #######################################################
	cout << "\n";		
	auto finish_total_time = std::chrono::high_resolution_clock::now();	
	std::chrono::duration<double> elapsed_total_time = finish_total_time - start_total_time;
	cout << "->-> Total Elapsed Time: " << elapsed_total_time.count() << " s\n";	
	cout << " \n";

	cout << "################################################################################" << "\n";
	cout << "## Stop - NuMagSANS ############################################################" << "\n";
	cout << "################################################################################" << "\n\n\n\n";

}
