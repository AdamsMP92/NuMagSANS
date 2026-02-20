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
#include "NuMagSANSlib_MagData.h"
#include "NuMagSANSlib_NucData.h"
#include "NuMagSANSlib_StructureData.h"
#include "NuMagSANSlib_SANSData.h"
#include "NuMagSANSlib_SpectralData.h"
#include "NuMagSANSlib_gpuKernel.h"

using namespace std;

void NuMagSANS_Calculator(InputFileData* InputData, \
					      NucDataProperties* NucDataProp,\
                          MagDataProperties* MagDataProp,\
                          StructDataProperties* StructDataProp, \
                          int Data_File_Index){

	cudaError_t err;
	
	LogSystem::write("################################################################################");
	LogSystem::write("## Run - NuMagSANS #############################################################");
	LogSystem::write("################################################################################");
	LogSystem::write("");

	// start time measurement #################################################################
	auto start_total_time = std::chrono::high_resolution_clock::now();

	// initialize nuclear data ################################################################
	NuclearData NucData, NucData_gpu;
	if(InputData->NucData_activate_flag){
		init_NuclearData(&NucData, &NucData_gpu, NucDataProp, InputData, Data_File_Index);
		cudaDeviceSynchronize();
		err = cudaGetLastError();
		if (err != cudaSuccess) {
			LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
		}
	}
	
	// initialize magnetization data ##########################################################
	MagnetizationData MagData, MagData_gpu;
	if(InputData->MagData_activate_flag){
		init_MagnetizationData(&MagData, &MagData_gpu, MagDataProp, InputData, Data_File_Index);
		cudaDeviceSynchronize();
		err = cudaGetLastError();
		if (err != cudaSuccess) {
			LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
		}
		//disp_MagnetizationData(&MagData);
	}

	// initialize structure data ##############################################################
	StructureData StructData, StructData_gpu;
	if(InputData->StructData_activate_flag){
		init_StructureData(&StructData, &StructData_gpu, StructDataProp, InputData);
		cudaDeviceSynchronize();
		err = cudaGetLastError();
		if (err != cudaSuccess) {
			LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
		}
		//disp_StructureData(&StructData);
	}

	// initialize scattering data #############################################################
	ScatteringData SANSData, SANSData_gpu;
	init_ScatteringData(InputData, &SANSData, &SANSData_gpu);
	cudaDeviceSynchronize();
	err = cudaGetLastError();
	if (err != cudaSuccess) {
		LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
	}

	// initialize spectral data ##############################################################
	SpectralData SpecData, SpecData_gpu;
	init_SpectralData(InputData, &SANSData, &SpecData, &SpecData_gpu);
	cudaDeviceSynchronize();
	err = cudaGetLastError();
	if (err != cudaSuccess) {
		LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
	}
	
	// initialize scaling factors #############################################################
	ScalingFactors ScalFactors;
	init_ScalingFactors(&ScalFactors, InputData, &MagData, &NucData, &SANSData);
	
	// compute 2D SANS cross sections #########################################################
	int L = (*SANSData.N_q) * (*SANSData.N_theta);
	LogSystem::write("total number of Fourier space bins: " + std::to_string(L));
	
		// Pure Magnetic Scattering Calculator without structure data #####################################################################
		if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 0 && InputData->StructData_activate_flag == 0){
			LogSystem::write("run: Atomistic_MagSANS_Kernel_dilute");
			Atomistic_MagSANS_Kernel_dilute<<<(L+255)/256, 256>>>(MagData_gpu, SANSData_gpu);
			cudaDeviceSynchronize();
			err = cudaGetLastError();
			if (err != cudaSuccess) {
				LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
			}
		}

		// Pure Nuclear Scattering Calculator without structure data #######################################################################
		if(InputData->MagData_activate_flag == 0 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 0){
			LogSystem::write("run: Atomistic_NucSANS_Kernel_dilute");
			Atomistic_NucSANS_Kernel_dilute<<<(L+255)/256, 256>>>(NucData_gpu, SANSData_gpu);
			cudaDeviceSynchronize();
			err = cudaGetLastError();
			if (err != cudaSuccess) {
				LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
			}
		}

		// Combined Magnetic and Nuclear Scattering Calculator without structure data ######################################################
		if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 0){
			LogSystem::write("run: Atomistic_NuMagSANS_Kernel_dilute");
			Atomistic_NuMagSANS_Kernel_dilute<<<(L+255)/256, 256>>>(NucData_gpu, MagData_gpu, SANSData_gpu);
			cudaDeviceSynchronize();
			err = cudaGetLastError();
			if (err != cudaSuccess) {
				LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
			}
		}


		// Pure Magnetic Scattering Calculator with structure data #########################################################################
		if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 0 && InputData->StructData_activate_flag == 1){
			LogSystem::write("run: Atomistic_MagSANS_Kernel");
			Atomistic_MagSANS_Kernel<<<(L+255)/256, 256>>>(MagData_gpu, StructData_gpu, SANSData_gpu);
			cudaDeviceSynchronize();
			err = cudaGetLastError();
			if (err != cudaSuccess) {
				LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
			}
		}

		// Pure Nuclear Scattering Calculator with structure data ###########################################################################
		if(InputData->MagData_activate_flag == 0 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 1){
			LogSystem::write("run: Atomistic_NucSANS_Kernel");
			Atomistic_NucSANS_Kernel<<<(L+255)/256, 256>>>(NucData_gpu, StructData_gpu, SANSData_gpu);
			cudaDeviceSynchronize();
			err = cudaGetLastError();
			if (err != cudaSuccess) {
				LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
			}
		}

		// Combined Magnetic and Nuclear Scattering Calculator with structure data #########################################################
		if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 1){
			LogSystem::write("run: Atomistic_NuMagSANS_Kernel");
			Atomistic_NuMagSANS_Kernel<<<(L+255)/256, 256>>>(NucData_gpu, MagData_gpu, StructData_gpu, SANSData_gpu);
			cudaDeviceSynchronize();
			err = cudaGetLastError();
			if (err != cudaSuccess) {
				LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
			}
		}

	// compute azimuthal average 1D ###########################################################
	bool compute_1D_azimuthal_average =
	       InputData->output_unpolarized_nuclear_SANS_cross_section_1D_flag
    		|| InputData->output_unpolarized_magnetic_SANS_cross_section_1D_flag
    		|| InputData->output_polarized_magnetic_SANS_cross_section_1D_flag
    		|| InputData->output_nuclear_magnetic_SANS_cross_section_1D_flag
    		|| InputData->output_spin_flip_magnetic_SANS_cross_section_1D_flag
    		|| InputData->output_chiral_magnetic_SANS_cross_section_1D_flag
    		|| InputData->output_pm_spin_flip_SANS_cross_section_1D_flag
    		|| InputData->output_mp_spin_flip_SANS_cross_section_1D_flag
    		|| InputData->output_pp_non_spin_flip_SANS_cross_section_1D_flag
    		|| InputData->output_mm_non_spin_flip_SANS_cross_section_1D_flag
    		|| InputData->output_m_sanspol_cross_section_1D_flag
    		|| InputData->output_p_sanspol_cross_section_1D_flag;

	if(compute_1D_azimuthal_average){

		LogSystem::write("run: azimuthal averaging");
		AzimuthalAverage<<<(L+255)/256, 256>>>(SANSData_gpu);
		cudaDeviceSynchronize();
		err = cudaGetLastError();
		if (err != cudaSuccess) {
		   LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err)); 
		}

	}

	// compute 1D pair distance distribution and correlation function #########################
	bool compute_1D_realspace =
	       InputData->output_nuclear_pair_distance_distribution_1D_flag
    		|| InputData->output_unpolarized_pair_distance_distribution_1D_flag
    		|| InputData->output_polarized_pair_distance_distribution_1D_flag
    		|| InputData->output_nuclear_magnetic_pair_distance_distribution_1D_flag
    		|| InputData->output_spin_flip_pair_distance_distribution_1D_flag
    		|| InputData->output_chiral_pair_distance_distribution_1D_flag
    		|| InputData->output_pm_spin_flip_pair_distance_distribution_1D_flag
    		|| InputData->output_mp_spin_flip_pair_distance_distribution_1D_flag
    		|| InputData->output_pp_non_spin_flip_pair_distance_distribution_1D_flag
    		|| InputData->output_mm_non_spin_flip_pair_distance_distribution_1D_flag
    		|| InputData->output_p_sanspol_pair_distance_distribution_1D_flag
    		|| InputData->output_m_sanspol_pair_distance_distribution_1D_flag
    		|| InputData->output_nuclear_correlation_function_1D_flag
    		|| InputData->output_unpolarized_correlation_function_1D_flag
    		|| InputData->output_polarized_correlation_function_1D_flag
    		|| InputData->output_nuclear_magnetic_correlation_function_1D_flag
    		|| InputData->output_spin_flip_correlation_function_1D_flag
    		|| InputData->output_chiral_correlation_function_1D_flag
    		|| InputData->output_pm_spin_flip_correlation_function_1D_flag
    		|| InputData->output_mp_spin_flip_correlation_function_1D_flag
    		|| InputData->output_pp_non_spin_flip_correlation_function_1D_flag
    		|| InputData->output_mm_non_spin_flip_correlation_function_1D_flag
    		|| InputData->output_p_sanspol_correlation_function_1D_flag
    		|| InputData->output_m_sanspol_correlation_function_1D_flag;

	if(compute_1D_realspace){
		LogSystem::write("run: 1D correlation functions");
		DistributionFunctions<<<(L+255), 256>>>(SANSData_gpu);
		cudaDeviceSynchronize();
		err = cudaGetLastError();
		if (err != cudaSuccess) {
	   		LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
		}
	}

	// compute 2D correlation functions ######################################################
	bool compute_2D_correlation =
	       InputData->output_nuclear_correlation_function_2D_flag
    		|| InputData->output_unpolarized_correlation_function_2D_flag
    		|| InputData->output_polarized_correlation_function_2D_flag
    		|| InputData->output_nuclear_magnetic_correlation_function_2D_flag
    		|| InputData->output_spin_flip_correlation_function_2D_flag
    		|| InputData->output_chiral_correlation_function_2D_flag
    		|| InputData->output_pm_spin_flip_correlation_function_2D_flag
    		|| InputData->output_mp_spin_flip_correlation_function_2D_flag
    		|| InputData->output_pp_non_spin_flip_correlation_function_2D_flag
    		|| InputData->output_mm_non_spin_flip_correlation_function_2D_flag
    		|| InputData->output_m_sanspol_correlation_function_2D_flag
    		|| InputData->output_p_sanspol_correlation_function_2D_flag;	


	if(compute_2D_correlation){
		LogSystem::write("run: 2D correlation functions");
		CorrelationFunction_2D<<<(L+255)/256, 256>>>(SANSData_gpu);
		cudaDeviceSynchronize();
		err = cudaGetLastError();
		if (err != cudaSuccess) {
		   LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
		}
	}


	// compute angular spectral intensities ###################################################
	if(InputData->AngularSpec_activate_flag){
	
		LogSystem::write("run: angular spectrum analyzer");
		ComputeSpectralDecomposition<<<(L+255)/256, 256>>>(SANSData_gpu, SpecData_gpu);
		err = cudaGetLastError();
		if (err != cudaSuccess) {
		   LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
		}

	// compute angular spectral amplitudes ####################################################
        	LogSystem::write("run: angular amplitude spectrum analyzer");
		ComputeAngularSpectrumAmplitudes<<<(L+255)/256, 256>>>(SpecData_gpu);
		err = cudaGetLastError();
        	if (err != cudaSuccess) {
           		LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
        	}

	}


	// copy scattering data from GPU to RAM ###################################################
	copyGPU2RAM_ScatteringData(&SANSData, &SANSData_gpu);

	// scaling of the scattering data on RAM ##################################################
	LogSystem::write("scaling of SANSdata...");
	scale_ScatteringData(&ScalFactors, &SANSData, InputData);

	// write scattering data to csv files #####################################################
	write2CSVtable_ScatteringData(InputData, &SANSData, Data_File_Index);

	if(InputData->AngularSpec_activate_flag){
		// copy spectral data from GPU to RAM #####################################################
		copyGPU2RAM_SpectralData(&SpecData, &SpecData_gpu);

		// scaling of the spectral data on RAM ####################################################
		scale_SpectralData(&ScalFactors, &SpecData, InputData);
	
		// write spectral data to csv files #######################################################
		write2CSV_SpectralData(InputData, &SpecData, Data_File_Index);
	}

	// free memory ############################################################################
	LogSystem::write("free memory...");
	if(InputData->NucData_activate_flag){
		free_NuclearData(&NucData, &NucData_gpu);
	}

	if(InputData->MagData_activate_flag){
		free_MagnetizationData(&MagData, &MagData_gpu);
	}

	if(InputData->StructData_activate_flag){
		free_StructureData(&StructData, &StructData_gpu);
	}
	
	// free scattering data
	free_ScatteringData(&SANSData, &SANSData_gpu);

	// free spectral data
	free_SpectralData(&SpecData, &SpecData_gpu);
	
	// print result of time measurement #######################################################
	LogSystem::write("");
	auto finish_total_time = std::chrono::high_resolution_clock::now();	
	std::chrono::duration<double> elapsed_total_time = finish_total_time - start_total_time;
	LogSystem::write("->-> Total Elapsed Time: " + std::to_string(elapsed_total_time.count()) + " s");
	LogSystem::write("");

	LogSystem::write("################################################################################");
	LogSystem::write("## Finished - NuMagSANS ########################################################");
	LogSystem::write("################################################################################");
	LogSystem::write("");
	LogSystem::write("");
}
