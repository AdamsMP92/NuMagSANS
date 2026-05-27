
.. _program_listing_file_src_NuMagSANSlib_ExportData.h:

Program Listing for File NuMagSANSlib_ExportData.h
==================================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_NuMagSANSlib_ExportData.h>` (``src/NuMagSANSlib_ExportData.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   
   // inline void ExportData(InputFileData* InputData,
   //                     ScalingFactors* ScalFactors,
   //                     NuMagSANSData* Data,
   //                     int Data_File_Index){
   
   //  // copy scattering data from GPU to RAM ###################################################
   //  copyGPU2RAM_ScatteringData(&Data->SANSData, &Data->SANSData_gpu);
   
   //  // scaling of the scattering data on RAM ##################################################
   //  LogSystem::write("scaling of SANSdata...");
   //  scale_ScatteringData(ScalFactors, &Data->SANSData, InputData);
   
   //  // write scattering data to csv files #####################################################
   //  write2CSVtable_ScatteringData(InputData, &Data->SANSData, Data_File_Index);
   
   //  if(InputData->AngularSpec_activate_flag){
   //      // copy spectral data from GPU to RAM #################################################
   //      copyGPU2RAM_SpectralData(&Data->SpecData, &Data->SpecData_gpu);
   
   //      // scaling of the spectral data on RAM ################################################
   //      scale_SpectralData(ScalFactors, &Data->SpecData, InputData);
       
   //      // write spectral data to csv files ###################################################
   //      write2CSV_SpectralData(InputData, &Data->SpecData, Data_File_Index);
   //  }
   // }
   
   inline void ExportData(InputFileData* InputData,
                          ScalingFactors* ScalFactors,
                          ScatteringData* SANSData, ScatteringData* SANSData_gpu,
                          SpectralData* SpecData, SpectralData* SpecData_gpu,
                          int Data_File_Index){
   
       // copy scattering data from GPU to RAM ###################################################
       copyGPU2RAM_ScatteringData(SANSData, SANSData_gpu);
   
       // scaling of the scattering data on RAM ##################################################
       LogSystem::write("scaling of SANSdata...");
       scale_ScatteringData(ScalFactors, SANSData, InputData);
   
       // write scattering data to csv files #####################################################
       write2CSVtable_ScatteringData(InputData, SANSData, Data_File_Index);
   
       if(InputData->AngularSpec_activate_flag){
           // copy spectral data from GPU to RAM #################################################
           copyGPU2RAM_SpectralData(SpecData, SpecData_gpu);
   
           // scaling of the spectral data on RAM ################################################
           scale_SpectralData(ScalFactors, SpecData, InputData);
       
           // write spectral data to csv files ###################################################
           write2CSV_SpectralData(InputData, SpecData, Data_File_Index);
       }
   }
