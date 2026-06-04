
.. _program_listing_file_src_wrapper_NuMagSANSlib_Data.h:

Program Listing for File NuMagSANSlib_Data.h
============================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_wrapper_NuMagSANSlib_Data.h>` (``src/wrapper/NuMagSANSlib_Data.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   
   struct NuMagSANSData {
       ScalingFactors ScalFactors;
       NuclearData NucData, NucData_gpu;
       MagnetizationData MagData, MagData_gpu;
       StructureData StructData, StructData_gpu;
       RotationData RotData, RotData_gpu;
       ScatteringData SANSData, SANSData_gpu;
       SpectralData SpecData, SpecData_gpu;
   };
