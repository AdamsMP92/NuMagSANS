
.. _program_listing_file_src_NuMagSANSlib_StringCompare.h:

Program Listing for File NuMagSANSlib_StringCompare.h
=====================================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_NuMagSANSlib_StringCompare.h>` (``src/NuMagSANSlib_StringCompare.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   // File         : NuMagSANSlib_StringCompare.h
   // Author       : Michael Philipp ADAMS, M.Sc.
   // Company      : University of Luxembourg
   // Department   : Department of Physics and Materials Sciences
   // Group        : NanoMagnetism Group
   // Group Leader : Prof. Andreas Michels
   // Version      : 14 November 2024
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
   
   void sort_integ(int *x, int K){
       int buff = 0;
       for(int i = 0; i < K; i++){
           for(int j = i + 1; j < K; j++){
               if(x[i] > x[j]){
                   buff = x[i];
                   x[i] = x[j];
                   x[j] = buff;
               }
           }
       }
   }
   
   
   bool Check_SingleUnderscore(std::string * Names, int Number_Of_Names){
   
       bool CheckFlag = true;
       int underscore_count = 0;
       string name;
       for(int k=0; k<Number_Of_Names; k++){
           underscore_count = 0;
           for (char ch : Names[k]) {
               if (ch == '_') {
                   underscore_count++;
               }
           }
           if(underscore_count > 1){
               CheckFlag = false;
               break;
           }
       }
       return CheckFlag;
   }
   
   bool Check_SinglePoint(std::string * Names, int Number_Of_Names){
   
       bool CheckFlag = true;
       int underscore_count = 0;
       string name;
       for(int k=0; k<Number_Of_Names; k++){
           underscore_count = 0;
           for (char ch : Names[k]) {
               if (ch == '.') {
                   underscore_count++;
               }
           }
           if(underscore_count > 1){
               CheckFlag = false;
               break;
           }
       }
       return CheckFlag;
   }
   
   
   
   void Split_Nom_Idx(std::string *Names, int Number_Of_Names, std::string*Nom, int*Idx){
   
       string name;
       for(int k = 0; k < Number_Of_Names; k++){
           name = Names[k];
   
           size_t underscore_pos = name.find('_');
   
           if(underscore_pos != std::string::npos){
               Nom[k] = name.substr(0, underscore_pos);
               Idx[k] = std::stoi(name.substr(underscore_pos + 1));
           }
           else{
               Nom[k] = name;
               Idx[k] = -1;
           }
       }
   }
   
   void Split_Nom_Idx_Type(std::string *Names, int Number_Of_Names, std::string*Nom, int*Idx, std::string*Type){
   
       string name;
       for(int k = 0; k < Number_Of_Names; k++){
           name = Names[k];
   
           size_t underscore_pos = name.find('_');
           size_t point_pos = name.find('.');
   
           if(underscore_pos != std::string::npos){
               Nom[k] = name.substr(0, underscore_pos);
               Idx[k] = std::stoi(name.substr(underscore_pos + 1, point_pos - underscore_pos - 1));
               Type[k] = name.substr(point_pos + 1);
           }
           else{
               Nom[k] = name;
               Idx[k] = -1;
               Type[k] = "";
           }
       }
   }
   
   bool Check_Nom_Consistency(std::string *Nom, int Number_Of_Names){
   
       std::string ReferenceNom = Nom[0];
       for(int i = 0; i < Number_Of_Names; i++){
           if(Nom[i] != ReferenceNom){
               return false;   
           }
       }
       return true;
   }
   
   bool Check_Idx_Consistency(int*Idx, int Number_Of_Names){
   
       for(int k = 0; k < Number_Of_Names; k++){
           if(Idx[k] != k+1){
               return false;
           }
       }
       return true;
   
   }
   
   
   void CheckStrings_FolderNames(std::string * Names, int Number_Of_Names, bool*CheckFlag, std::string*Nom){
   
       bool UnderscoreCheckFlag = Check_SingleUnderscore(Names, Number_Of_Names);
   
       if(UnderscoreCheckFlag){
           string Nom_array[Number_Of_Names];
           int Idx_array[Number_Of_Names];
           Split_Nom_Idx(Names, Number_Of_Names, Nom_array, Idx_array);            
           sort_integ(Idx_array, Number_Of_Names);
   
           bool NomConsistencyFlag = Check_Nom_Consistency(Nom_array, Number_Of_Names);        
           bool IdxConsistencyFlag = Check_Idx_Consistency(Idx_array, Number_Of_Names);
           if(NomConsistencyFlag && IdxConsistencyFlag){
               *Nom = Nom_array[0];
               *CheckFlag = true;
           }
           else{
               *Nom = "";
               *CheckFlag = false;
           }
       }
       else{
                   *Nom = "";
                   *CheckFlag = false;
       }
   }
   
   
   void CheckStrings_ElementNames(std::string * Names, int Number_Of_Names, bool*CheckFlag, std::string*Nom, std::string*Type){
   
       bool UnderscoreCheckFlag = Check_SingleUnderscore(Names, Number_Of_Names);
       bool PointCheckFlag = Check_SinglePoint(Names, Number_Of_Names);
   
       if(UnderscoreCheckFlag && PointCheckFlag){
           std::string Nom_array[Number_Of_Names];
           int Idx_array[Number_Of_Names];
           std::string Type_array[Number_Of_Names];
           Split_Nom_Idx_Type(Names, Number_Of_Names, Nom_array, Idx_array, Type_array);
           sort_integ(Idx_array, Number_Of_Names);
   
           bool NomConsistencyFlag = Check_Nom_Consistency(Nom_array, Number_Of_Names);        
           bool IdxConsistencyFlag = Check_Idx_Consistency(Idx_array, Number_Of_Names);
           bool TypeConsistencyFlag = Check_Nom_Consistency(Type_array, Number_Of_Names);      
   
           if(NomConsistencyFlag && IdxConsistencyFlag && TypeConsistencyFlag){
               *Nom = Nom_array[0];
               *Type = Type_array[0];
               *CheckFlag = true;
           }
           else{
               *Nom = "";
               *Type = "";
               *CheckFlag = false;
           }
       }
       else{
           *Nom = "";
           *Type = "";
           *CheckFlag = false;
       }   
   }
   
   
