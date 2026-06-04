
.. _program_listing_file_src_helper_NuMagSANSlib_TimeMeasure.h:

Program Listing for File NuMagSANSlib_TimeMeasure.h
===================================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_helper_NuMagSANSlib_TimeMeasure.h>` (``src/helper/NuMagSANSlib_TimeMeasure.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   
   struct TimeMeasure {
       std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
   };
   
   inline TimeMeasure StartTimeMeasure(){
       TimeMeasure Time = {};
       Time.start_time = std::chrono::high_resolution_clock::now();
       return Time;
   }
   
   inline void LogElapsedTime(const TimeMeasure& Time){
       auto finish_time = std::chrono::high_resolution_clock::now();
       std::chrono::duration<double> elapsed_time = finish_time - Time.start_time;
       LogSystem::write("");
       LogSystem::write("->-> Total Elapsed Time: " + std::to_string(elapsed_time.count()) + " s");
       LogSystem::write("");
   }
