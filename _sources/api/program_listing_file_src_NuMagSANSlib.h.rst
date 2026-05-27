
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
   #include "NuMagSANSlib_Data.h"
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
   
       NuMagSANSData Data;
   
           // ScalingFactors ScalFactors;
           // NuclearData NucData, NucData_gpu;
           // MagnetizationData MagData, MagData_gpu;
           // StructureData StructData, StructData_gpu;
           // RotationData RotData, RotData_gpu;
           // ScatteringData SANSData, SANSData_gpu;
           // SpectralData SpecData, SpecData_gpu;
   
       // InitializeData(InputData,
       //             NucDataProp,
       //             MagDataProp,
       //             StructDataProp,
       //             RotDataProp,
       //             &NucData, &NucData_gpu,
       //             &MagData, &MagData_gpu,
       //             &StructData, &StructData_gpu,
       //             &RotData, &RotData_gpu,
       //             &SANSData, &SANSData_gpu,
       //             &SpecData, &SpecData_gpu,
       //             Data_File_Index);
       InitializeData(InputData,
                      NucDataProp,
                      MagDataProp,
                      StructDataProp,
                      RotDataProp,
                      &Data, 
                      Data_File_Index);
   
       LogCurrentGPUMemoryDifference(MemoryBeforeRun);
       
           // initialize scaling factors #############################################################
           // init_ScalingFactors(&ScalFactors, InputData, &MagData, &NucData, &SANSData);
       
           // SelectKernelRun(InputData,
           //              &Data.NucData_gpu,
           //              &Data.MagData_gpu,
           //              &Data.StructData_gpu,
           //              &Data.RotData_gpu,
           //              &Data.SANSData,
           //              &Data.SANSData_gpu);
           
           SelectKernelRun(InputData, &Data);
   
           // KernelPostprocessRun(InputData,
           //                   &SANSData,
           //                   &SANSData_gpu,
           //                   &SpecData_gpu);
   
           KernelPostprocessRun(InputData, &Data);
   
       // copy scattering data from GPU to RAM ###################################################
       // ExportData(InputData,
       //         &ScalFactors,
       //         &SANSData, &SANSData_gpu,
       //         &SpecData, &SpecData_gpu,
       //         Data_File_Index);
       ExportData(InputData, &Data, Data_File_Index);
   
       // free memory ############################################################################
       // FreeData(InputData,
       //       &MagData, &MagData_gpu,
       //       &NucData, &NucData_gpu,
       //       &StructData, &StructData_gpu,
       //       &RotData, &RotData_gpu,
       //       &SANSData, &SANSData_gpu,
       //       &SpecData, &SpecData_gpu);
       FreeData(InputData, &Data);
       
       // print result of time measurement #######################################################
       LogElapsedTime(TotalTime);
   
       LogCurrentGPUMemoryDifference(MemoryBeforeRun);
   
       LogSystem::write("################################################################################");
       LogSystem::write("## Finished - NuMagSANS ########################################################");
       LogSystem::write("################################################################################");
       LogSystem::write("");
       LogSystem::write("");
   }
