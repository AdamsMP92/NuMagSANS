
.. _program_listing_file_src_helper_NuMagSANSlib_MemoryInfo.h:

Program Listing for File NuMagSANSlib_MemoryInfo.h
==================================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_helper_NuMagSANSlib_MemoryInfo.h>` (``src/helper/NuMagSANSlib_MemoryInfo.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   
   struct GPUMemoryInfo {
       size_t free_bytes;
       size_t total_bytes;
       double free_mb;
       double used_mb;
   };
   
   inline GPUMemoryInfo GetGPUMemoryInfo() {
       GPUMemoryInfo Info = {};
       cudaMemGetInfo(&Info.free_bytes, &Info.total_bytes);
       Info.free_mb = Info.free_bytes / 1024.0 / 1024.0;
       Info.used_mb = (Info.total_bytes - Info.free_bytes) / 1024.0 / 1024.0;
       return Info;
   }
   
   inline void LogCurrentGPUMemoryDifference(const GPUMemoryInfo& MemoryReference) {
       GPUMemoryInfo MemoryCurrent = GetGPUMemoryInfo();
       double allocated_mb = MemoryCurrent.used_mb - MemoryReference.used_mb;
       LogSystem::write("");
       LogSystem::write("GPU Memory Check: allocated bytes: " + std::to_string(allocated_mb) +
                        " MB, used bytes: " + std::to_string(MemoryCurrent.used_mb) +
                        " MB, free bytes: " + std::to_string(MemoryCurrent.free_mb) + " MB");
       LogSystem::write("");
   }
