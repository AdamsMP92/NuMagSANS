
.. _program_listing_file_src_NuMagSANSlib_CUDAError.h:

Program Listing for File NuMagSANSlib_CUDAError.h
=================================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_NuMagSANSlib_CUDAError.h>` (``src/NuMagSANSlib_CUDAError.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   
   inline bool CheckCUDA(cudaError_t err,
                         const std::string& context){
       if(err != cudaSuccess){
           LogSystem::write("CUDA error in " + context + ": " + cudaGetErrorString(err));
           return false;
       }
   
       return true;
   }
   
   inline bool CheckCUDALastError(const std::string& context){
       bool success = true;
   
       success = CheckCUDA(cudaGetLastError(),
                           context) && success;
   
       success = CheckCUDA(cudaDeviceSynchronize(),
                           context + " synchronization") && success;
   
       return success;
   }
   
   inline bool CheckCUDAKernelRun(const std::string& kernel_name){
       bool success = true;
   
       success = CheckCUDA(cudaGetLastError(),
                           "kernel launch " + kernel_name) && success;
   
       success = CheckCUDA(cudaDeviceSynchronize(),
                           "kernel execution " + kernel_name) && success;
   
       return success;
   }
