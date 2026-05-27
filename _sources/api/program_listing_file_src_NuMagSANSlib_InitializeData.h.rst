
.. _program_listing_file_src_NuMagSANSlib_InitializeData.h:

Program Listing for File NuMagSANSlib_InitializeData.h
======================================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_NuMagSANSlib_InitializeData.h>` (``src/NuMagSANSlib_InitializeData.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   
   inline void InitializeData(InputFileData* InputData,
                              NucDataProperties* NucDataProp,
                              MagDataProperties* MagDataProp,
                              StructDataProperties* StructDataProp,
                              RotDataProperties* RotDataProp,
                              NuMagSANSData* Data,
                              int Data_File_Index){
   
       cudaError_t err;
   
       // initialize nuclear data ################################################################
       if(InputData->NucData_activate_flag){
           init_NuclearData(&Data->NucData, &Data->NucData_gpu, NucDataProp, InputData, Data_File_Index);
           cudaDeviceSynchronize();
           err = cudaGetLastError();
           CheckKernelLaunch(err);
       }
       
       // initialize magnetization data ##########################################################
       if(InputData->MagData_activate_flag){
           init_MagnetizationData(&Data->MagData, &Data->MagData_gpu, MagDataProp, InputData, Data_File_Index);
           cudaDeviceSynchronize();
           err = cudaGetLastError();
           CheckKernelLaunch(err);
           //disp_MagnetizationData(&Data->MagData);
       }
   
       // initialize structure data ##############################################################
       if(InputData->StructData_activate_flag){
           init_StructureData(&Data->StructData, &Data->StructData_gpu, StructDataProp, InputData);
           cudaDeviceSynchronize();
           err = cudaGetLastError();
           CheckKernelLaunch(err);
           //disp_StructureData(&Data->StructData);
       }
   
       // initialize rotation data ###############################################################
       if(InputData->RotData_activate_flag){
           init_RotationData(&Data->RotData, &Data->RotData_gpu, RotDataProp, InputData);
           cudaDeviceSynchronize();
           err = cudaGetLastError();
           CheckKernelLaunch(err);
       }
   
       // initialize scattering data #############################################################
       init_ScatteringData(InputData, &Data->SANSData, &Data->SANSData_gpu);
       cudaDeviceSynchronize();
       err = cudaGetLastError();
       CheckKernelLaunch(err);
   
       // initialize spectral data ###############################################################
       init_SpectralData(InputData, &Data->SANSData, &Data->SpecData, &Data->SpecData_gpu);
       cudaDeviceSynchronize();
       err = cudaGetLastError();
       CheckKernelLaunch(err);
   }
