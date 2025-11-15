
.. _program_listing_file_src_NuMagSANSlib_HelperFun.h:

Program Listing for File NuMagSANSlib_HelperFun.h
=================================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_NuMagSANSlib_HelperFun.h>` (``src/NuMagSANSlib_HelperFun.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   // File         : NuMagSANSlib_ReadWrite.h
   // Author       : Michael Philipp ADAMS, M.Sc.
   // Company      : University of Luxembourg
   // Department   : Department of Physics and Materials Sciences
   // Group        : NanoMagnetism Group
   // Group Leader : Prof. Andreas Michels
   // Version      : 22 November 2024
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
   
   using namespace std;
   
   
   
   
   unsigned long int SumIntList2(int** List, int NumberOfColumns){
   
       unsigned long int L = 0;
   
       for(int k = 0; k < NumberOfColumns; k++){
           L += *List[k];
       }
   
       return L;
   
   }
   
