
.. _program_listing_file_src_NuMagSANSlib_StructureDataExplorer.h:

Program Listing for File NuMagSANSlib_StructureDataExplorer.h
=============================================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_NuMagSANSlib_StructureDataExplorer.h>` (``src/NuMagSANSlib_StructureDataExplorer.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   // File         : NuMagSANSlib_MagDataObserver.h
   // Author       : Michael Philipp ADAMS, M.Sc.
   // Company      : University of Luxembourg
   // Department   : Department of Physics and Materials Sciences
   // Group        : NanoMagnetism Group
   // Group Leader : Prof. Andreas Michels
   // Version      : 23 October 2024
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
   #include <limits.h>
   
   using namespace std;
   
   struct StructDataProperties{
   
       string GlobalFilePath;
   
       int Number_Of_Elements;
       
   };
   
   // ###############################################################################################################################################
   // helper functions ##############################################################################################################################
   // ###############################################################################################################################################
   
   void get_GlobalStructDataPath(std::string Local_StructDataPath, StructDataProperties* StructDataProp){
   
       char tmp[PATH_MAX];
       getcwd(tmp, PATH_MAX);  // Get the current working directory
       std::string tmp_string = tmp;
       StructDataProp->GlobalFilePath = tmp_string + "/" + Local_StructDataPath;
       cout << "found Global StructDataPath: " << StructDataProp->GlobalFilePath << "\n\n";
       
   }
   
   void NumberOfEntriesInStructureFile(int *NumberOfColumns, string filename){
   
       ifstream fin;
       fin.open(filename);
       std::string line;
   
       float x, y, z;
   
       int line_counter = 0;
       int error_counter = 0;
   
       while(std::getline(fin, line)){
           std::istringstream ss(line);
           if(ss >> x >> y >> z){
               line_counter += 1;
           } else{
               error_counter ++;
           }
       }
       fin.close();
       *NumberOfColumns = line_counter;
   
   }
   
   // Routine that checks number of subfolders in MagData directory
   bool StructData_Observer(std::string Local_StructDataPath, StructDataProperties*StructDataProp){
   
       cout << "##########################################################################################" << "\n";
       cout << "## Run - Structure File Explorer #########################################################" << "\n";
       cout << "##########################################################################################" << "\n\n";
   
       bool CheckFlag = false;
   
       // get global path of the MagData directory
       get_GlobalStructDataPath(Local_StructDataPath, StructDataProp);
   
       // count number of entries in the structure data file
       NumberOfEntriesInStructureFile(&StructDataProp->Number_Of_Elements, StructDataProp->GlobalFilePath);
       cout << "Number of Entries: " << StructDataProp->Number_Of_Elements << "\n\n";
   
   
       cout << "##########################################################################################" << "\n";
       cout << "## Stop - Structure File Explorer ########################################################" << "\n";
       cout << "##########################################################################################" << "\n\n\n\n";
   
       if(StructDataProp->Number_Of_Elements != 0){
           CheckFlag = true;
       }
   
       return CheckFlag;
   
   }
   
