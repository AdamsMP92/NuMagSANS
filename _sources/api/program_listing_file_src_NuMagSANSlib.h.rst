
.. _program_listing_file_src_NuMagSANSlib.h:

Program Listing for File NuMagSANSlib.h
=======================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_NuMagSANSlib.h>` (``src/NuMagSANSlib.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

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
   #include "NuMagSANSlib_gpuKernel.h"
   
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
   
       // start time measurement #################################################################
       auto start_total_time = std::chrono::high_resolution_clock::now();
   
       // record GPU memory state before data initialization ######################################
       size_t free_bytes_before_data_load, total_bytes_before_data_load;
       size_t free_bytes_after_data_load, total_bytes_after_data_load;
       cudaMemGetInfo(&free_bytes_before_data_load, &total_bytes_before_data_load);
       double used_mb_before_data_load = (total_bytes_before_data_load - free_bytes_before_data_load) / 1024.0 / 1024.0;
   
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
   
       // initialize rotation data ##############################################################
       RotationData RotData, RotData_gpu;
       if(InputData->RotData_activate_flag){
           init_RotationData(&RotData, &RotData_gpu, RotDataProp, InputData);
           cudaDeviceSynchronize();
           err = cudaGetLastError();
           if (err != cudaSuccess) {
               LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
           }
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
   
       // report GPU memory state after data initialization #######################################
       cudaMemGetInfo(&free_bytes_after_data_load, &total_bytes_after_data_load);
       double used_mb_after_data_load = (total_bytes_after_data_load - free_bytes_after_data_load) / 1024.0 / 1024.0;
       double loaded_data_mb = used_mb_after_data_load - used_mb_before_data_load;
       LogSystem::write("GPU Memory Check after data load: loaded data: " + std::to_string(loaded_data_mb) + " MB, used bytes: " + std::to_string(used_mb_after_data_load) + " MB, free bytes: " + std::to_string(free_bytes_after_data_load / 1024.0 / 1024.0) + " MB");
       
       // initialize scaling factors #############################################################
       ScalingFactors ScalFactors;
       init_ScalingFactors(&ScalFactors, InputData, &MagData, &NucData, &SANSData);
       
       // compute 2D SANS cross sections #########################################################
       int L = (*SANSData.N_q) * (*SANSData.N_theta);
       LogSystem::write("total number of Fourier space bins: " + std::to_string(L));
       
               // Pure Magnetic Scattering Calculator without structure data #####################################################################
               if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 0 && InputData->StructData_activate_flag == 0 && InputData->RotData_activate_flag == 0){
               LogSystem::write("run: Atomistic_MagSANS_Kernel_dilute");
               Atomistic_MagSANS_Kernel_dilute<<<(L+255)/256, 256>>>(MagData_gpu, SANSData_gpu);
               cudaDeviceSynchronize();
               err = cudaGetLastError();
               if (err != cudaSuccess) {
                   LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
               }
           }
   
               // Pure Nuclear Scattering Calculator without structure data #######################################################################
               if(InputData->MagData_activate_flag == 0 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 0 && InputData->RotData_activate_flag == 0){
               LogSystem::write("run: Atomistic_NucSANS_Kernel_dilute");
               Atomistic_NucSANS_Kernel_dilute<<<(L+255)/256, 256>>>(NucData_gpu, SANSData_gpu);
               cudaDeviceSynchronize();
               err = cudaGetLastError();
               if (err != cudaSuccess) {
                   LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
               }
           }
   
               // Combined Magnetic and Nuclear Scattering Calculator without structure data ######################################################
               if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 0 && InputData->RotData_activate_flag == 0){
               LogSystem::write("run: Atomistic_NuMagSANS_Kernel_dilute");
               Atomistic_NuMagSANS_Kernel_dilute<<<(L+255)/256, 256>>>(NucData_gpu, MagData_gpu, SANSData_gpu);
               cudaDeviceSynchronize();
               err = cudaGetLastError();
               if (err != cudaSuccess) {
                   LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
               }
           }
   
   
               // Pure Magnetic Scattering Calculator with structure data #########################################################################
               if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 0 && InputData->StructData_activate_flag == 1 && InputData->RotData_activate_flag == 0){
               LogSystem::write("run: Atomistic_MagSANS_Kernel");
               Atomistic_MagSANS_Kernel<<<(L+255)/256, 256>>>(MagData_gpu, StructData_gpu, SANSData_gpu);
               cudaDeviceSynchronize();
               err = cudaGetLastError();
               if (err != cudaSuccess) {
                   LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
               }
           }
   
               // Pure Nuclear Scattering Calculator with structure data ###########################################################################
               if(InputData->MagData_activate_flag == 0 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 1 && InputData->RotData_activate_flag == 0){
               LogSystem::write("run: Atomistic_NucSANS_Kernel");
               Atomistic_NucSANS_Kernel<<<(L+255)/256, 256>>>(NucData_gpu, StructData_gpu, SANSData_gpu);
               cudaDeviceSynchronize();
               err = cudaGetLastError();
               if (err != cudaSuccess) {
                   LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
               }
           }
   
               // Combined Magnetic and Nuclear Scattering Calculator with structure data #########################################################
               if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 1 && InputData->RotData_activate_flag == 0){
               LogSystem::write("run: Atomistic_NuMagSANS_Kernel");
               Atomistic_NuMagSANS_Kernel<<<(L+255)/256, 256>>>(NucData_gpu, MagData_gpu, StructData_gpu, SANSData_gpu);
               cudaDeviceSynchronize();
               err = cudaGetLastError();
               if (err != cudaSuccess) {
                   LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
               }
               }
   
               // Pure Magnetic Scattering Calculator with rotation data and without structure data ##############################################
               if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 0 && InputData->StructData_activate_flag == 0 && InputData->RotData_activate_flag == 1){
                   LogSystem::write("run: Atomistic_MagSANS_Kernel_RotDilute");
                   Atomistic_MagSANS_Kernel_RotDilute<<<(L+255)/256, 256>>>(MagData_gpu, RotData_gpu, SANSData_gpu);
                   cudaDeviceSynchronize();
                   err = cudaGetLastError();
                   if (err != cudaSuccess) {
                       LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
                   }
               }
   
               // Pure Nuclear Scattering Calculator with rotation data and without structure data ###############################################
               if(InputData->MagData_activate_flag == 0 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 0 && InputData->RotData_activate_flag == 1){
                   LogSystem::write("run: Atomistic_NucSANS_Kernel_RotDilute");
                   Atomistic_NucSANS_Kernel_RotDilute<<<(L+255)/256, 256>>>(NucData_gpu, RotData_gpu, SANSData_gpu);
                   cudaDeviceSynchronize();
                   err = cudaGetLastError();
                   if (err != cudaSuccess) {
                       LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
                   }
               }
   
               // Combined Magnetic and Nuclear Scattering Calculator with rotation data and without structure data ###############################
               if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 0 && InputData->RotData_activate_flag == 1){
                   LogSystem::write("run: Atomistic_NuMagSANS_Kernel_RotDilute");
                   Atomistic_NuMagSANS_Kernel_RotDilute<<<(L+255)/256, 256>>>(NucData_gpu, MagData_gpu, RotData_gpu, SANSData_gpu);
                   cudaDeviceSynchronize();
                   err = cudaGetLastError();
                   if (err != cudaSuccess) {
                       LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
                   }
               }
   
               // Pure Magnetic Scattering Calculator with structure and rotation data ###########################################################
               if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 0 && InputData->StructData_activate_flag == 1 && InputData->RotData_activate_flag == 1){
                   LogSystem::write("run: Atomistic_MagSANS_Kernel_StructRot");
                   Atomistic_MagSANS_Kernel_StructRot<<<(L+255)/256, 256>>>(MagData_gpu, StructData_gpu, RotData_gpu, SANSData_gpu);
                   cudaDeviceSynchronize();
                   err = cudaGetLastError();
                   if (err != cudaSuccess) {
                       LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
                   }
               }
   
               // Pure Nuclear Scattering Calculator with structure and rotation data ############################################################
               if(InputData->MagData_activate_flag == 0 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 1 && InputData->RotData_activate_flag == 1){
                   LogSystem::write("run: Atomistic_NucSANS_Kernel_StructRot");
                   Atomistic_NucSANS_Kernel_StructRot<<<(L+255)/256, 256>>>(NucData_gpu, StructData_gpu, RotData_gpu, SANSData_gpu);
                   cudaDeviceSynchronize();
                   err = cudaGetLastError();
                   if (err != cudaSuccess) {
                       LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
                   }
               }
   
               // Combined Magnetic and Nuclear Scattering Calculator with structure and rotation data ############################################
               if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 1 && InputData->RotData_activate_flag == 1){
                   LogSystem::write("run: Atomistic_NuMagSANS_Kernel_StructRot");
                   Atomistic_NuMagSANS_Kernel_StructRot<<<(L+255)/256, 256>>>(NucData_gpu, MagData_gpu, StructData_gpu, RotData_gpu, SANSData_gpu);
                   cudaDeviceSynchronize();
                   err = cudaGetLastError();
                   if (err != cudaSuccess) {
                       LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
                   }
               }
   
           // compute azimuthal average 1D ###########################################################
       bool compute_1D_azimuthal_average = any_active(InputData->OutFlags.SANS1D);
       bool compute_1D_corr = any_active(InputData->OutFlags.Corr1D);
       bool compute_1D_pair = any_active(InputData->OutFlags.PairDist1D);
   
       if(compute_1D_azimuthal_average || compute_1D_corr || compute_1D_pair){
   
           LogSystem::write("run: azimuthal averaging");
           AzimuthalAverage<<<(L+255)/256, 256>>>(SANSData_gpu);
           cudaDeviceSynchronize();
           err = cudaGetLastError();
           if (err != cudaSuccess) {
              LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err)); 
           }
   
       }
   
       // compute 1D pair distance distribution and correlation function #########################
       if(compute_1D_corr || compute_1D_pair){
           LogSystem::write("run: 1D correlation functions");
           DistributionFunctions<<<(L+255), 256>>>(SANSData_gpu);
           cudaDeviceSynchronize();
           err = cudaGetLastError();
           if (err != cudaSuccess) {
               LogSystem::write(std::string("kernel launch failed: ") + cudaGetErrorString(err));
           }
       }
   
       // compute 2D correlation functions ######################################################
       bool compute_2D_correlation = any_active(InputData->OutFlags.Corr2D);   
   
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
   
       if(InputData->RotData_activate_flag){
           free_RotationData(&RotData, &RotData_gpu);
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
