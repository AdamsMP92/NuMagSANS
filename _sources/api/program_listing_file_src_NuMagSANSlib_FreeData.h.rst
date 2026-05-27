
.. _program_listing_file_src_NuMagSANSlib_FreeData.h:

Program Listing for File NuMagSANSlib_FreeData.h
================================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_NuMagSANSlib_FreeData.h>` (``src/NuMagSANSlib_FreeData.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   
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
