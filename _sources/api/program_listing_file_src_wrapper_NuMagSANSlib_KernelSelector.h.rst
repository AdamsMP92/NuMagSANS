
.. _program_listing_file_src_wrapper_NuMagSANSlib_KernelSelector.h:

Program Listing for File NuMagSANSlib_KernelSelector.h
======================================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_wrapper_NuMagSANSlib_KernelSelector.h>` (``src/wrapper/NuMagSANSlib_KernelSelector.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   
   inline void SelectKernelRunMulti(InputFileData* InputData,
                                   NuclearData* NucData_gpu,
                                   MagnetizationData* MagData_gpu,
                                   StructureData* StructData_gpu,
                                   RotationData* RotData_gpu,
                                   ScatteringData* SANSData,
                                   ScatteringData* SANSData_gpu){
   
       int L = (*SANSData->N_q) * (*SANSData->N_theta);
       LogSystem::write("total number of Fourier space bins: " + std::to_string(L));
   
       // Pure Magnetic Scattering Calculator without structure data #############################
       if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 0 && InputData->StructData_activate_flag == 0 && InputData->RotData_activate_flag == 0){
           LogSystem::write("run: Atomistic_MagSANS_Kernel_dilute");
           Atomistic_MagSANS_Kernel_dilute<<<(L+255)/256, 256>>>(*MagData_gpu, *SANSData_gpu);
           CheckCUDAKernelRun("Atomistic_MagSANS_Kernel_dilute");
       }
   
       // Pure Nuclear Scattering Calculator without structure data ##############################
       if(InputData->MagData_activate_flag == 0 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 0 && InputData->RotData_activate_flag == 0){
           LogSystem::write("run: Atomistic_NucSANS_Kernel_dilute");
           Atomistic_NucSANS_Kernel_dilute<<<(L+255)/256, 256>>>(*NucData_gpu, *SANSData_gpu);
           CheckCUDAKernelRun("Atomistic_NucSANS_Kernel_dilute");
       }
   
       // Combined Magnetic and Nuclear Scattering Calculator without structure data #############
       if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 0 && InputData->RotData_activate_flag == 0){
           LogSystem::write("run: Atomistic_NuMagSANS_Kernel_dilute");
           Atomistic_NuMagSANS_Kernel_dilute<<<(L+255)/256, 256>>>(*NucData_gpu, *MagData_gpu, *SANSData_gpu);
           CheckCUDAKernelRun("Atomistic_NuMagSANS_Kernel_dilute");
       }
   
       // Pure Magnetic Scattering Calculator with structure data ################################
       if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 0 && InputData->StructData_activate_flag == 1 && InputData->RotData_activate_flag == 0){
           LogSystem::write("run: Atomistic_MagSANS_Kernel");
           Atomistic_MagSANS_Kernel<<<(L+255)/256, 256>>>(*MagData_gpu, *StructData_gpu, *SANSData_gpu);
           CheckCUDAKernelRun("Atomistic_MagSANS_Kernel");
       }
   
       // Pure Nuclear Scattering Calculator with structure data #################################
       if(InputData->MagData_activate_flag == 0 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 1 && InputData->RotData_activate_flag == 0){
           LogSystem::write("run: Atomistic_NucSANS_Kernel");
           Atomistic_NucSANS_Kernel<<<(L+255)/256, 256>>>(*NucData_gpu, *StructData_gpu, *SANSData_gpu);
           CheckCUDAKernelRun("Atomistic_NucSANS_Kernel");
       }
   
       // Combined Magnetic and Nuclear Scattering Calculator with structure data ################
       if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 1 && InputData->RotData_activate_flag == 0){
           LogSystem::write("run: Atomistic_NuMagSANS_Kernel");
           Atomistic_NuMagSANS_Kernel<<<(L+255)/256, 256>>>(*NucData_gpu, *MagData_gpu, *StructData_gpu, *SANSData_gpu);
           CheckCUDAKernelRun("Atomistic_NuMagSANS_Kernel");
       }
   
       // Pure Magnetic Scattering Calculator with rotation data and without structure data ######
       if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 0 && InputData->StructData_activate_flag == 0 && InputData->RotData_activate_flag == 1){
           LogSystem::write("run: Atomistic_MagSANS_Kernel_RotDilute");
           Atomistic_MagSANS_Kernel_RotDilute<<<(L+255)/256, 256>>>(*MagData_gpu, *RotData_gpu, *SANSData_gpu);
           CheckCUDAKernelRun("Atomistic_MagSANS_Kernel_RotDilute");
       }
   
       // Pure Nuclear Scattering Calculator with rotation data and without structure data #######
       if(InputData->MagData_activate_flag == 0 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 0 && InputData->RotData_activate_flag == 1){
           LogSystem::write("run: Atomistic_NucSANS_Kernel_RotDilute");
           Atomistic_NucSANS_Kernel_RotDilute<<<(L+255)/256, 256>>>(*NucData_gpu, *RotData_gpu, *SANSData_gpu);
           CheckCUDAKernelRun("Atomistic_NucSANS_Kernel_RotDilute");
       }
   
       // Combined Magnetic and Nuclear Scattering Calculator with rotation data #################
       if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 0 && InputData->RotData_activate_flag == 1){
           LogSystem::write("run: Atomistic_NuMagSANS_Kernel_RotDilute");
           Atomistic_NuMagSANS_Kernel_RotDilute<<<(L+255)/256, 256>>>(*NucData_gpu, *MagData_gpu, *RotData_gpu, *SANSData_gpu);
           CheckCUDAKernelRun("Atomistic_NuMagSANS_Kernel_RotDilute");
       }
   
       // Pure Magnetic Scattering Calculator with structure and rotation data ###################
       if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 0 && InputData->StructData_activate_flag == 1 && InputData->RotData_activate_flag == 1){
           LogSystem::write("run: Atomistic_MagSANS_Kernel_StructRot");
           Atomistic_MagSANS_Kernel_StructRot<<<(L+255)/256, 256>>>(*MagData_gpu, *StructData_gpu, *RotData_gpu, *SANSData_gpu);
           CheckCUDAKernelRun("Atomistic_MagSANS_Kernel_StructRot");
       }
   
       // Pure Nuclear Scattering Calculator with structure and rotation data ####################
       if(InputData->MagData_activate_flag == 0 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 1 && InputData->RotData_activate_flag == 1){
           LogSystem::write("run: Atomistic_NucSANS_Kernel_StructRot");
           Atomistic_NucSANS_Kernel_StructRot<<<(L+255)/256, 256>>>(*NucData_gpu, *StructData_gpu, *RotData_gpu, *SANSData_gpu);
           CheckCUDAKernelRun("Atomistic_NucSANS_Kernel_StructRot");
       }
   
       // Combined Magnetic and Nuclear Scattering Calculator with structure and rotation data ###
       if(InputData->MagData_activate_flag == 1 && InputData->NucData_activate_flag == 1 && InputData->StructData_activate_flag == 1 && InputData->RotData_activate_flag == 1){
           LogSystem::write("run: Atomistic_NuMagSANS_Kernel_StructRot");
           Atomistic_NuMagSANS_Kernel_StructRot<<<(L+255)/256, 256>>>(*NucData_gpu, *MagData_gpu, *StructData_gpu, *RotData_gpu, *SANSData_gpu);
           CheckCUDAKernelRun("Atomistic_NuMagSANS_Kernel_StructRot");
       }
   }
   
   
   inline void SelectKernelRun(InputFileData* InputData, 
                               NuMagSANSData* Data){
   
   
       SelectKernelRunMulti(InputData,
                            &Data->NucData_gpu,
                            &Data->MagData_gpu,
                            &Data->StructData_gpu,
                            &Data->RotData_gpu,
                            &Data->SANSData,
                            &Data->SANSData_gpu);
   
   }
