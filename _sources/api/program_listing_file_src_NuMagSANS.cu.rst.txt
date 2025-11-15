
.. _program_listing_file_src_NuMagSANS.cu:

Program Listing for File NuMagSANS.cu
=====================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_NuMagSANS.cu>` (``src/NuMagSANS.cu``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   // File         : NuMagSANS.cu
   // Author       : Dr. Michael Philipp ADAMS
   // Company      : University of Luxembourg
   // Department   : Department of Physics and Materials Sciences
   // Group        : NanoMagnetism Group
   // Group Leader : Prof. Andreas Michels
   // Version      : 16 October 2025
   // OS           : Linux Ubuntu
   // Language     : CUDA C++
   
   #include "NuMagSANSlib.h"
   
   using namespace std;
   
   int main(int argc, char* argv[]){
       
       // Input File Name
       string InputFileName = argv[1];
   
       // Initialize LogFileSystem
       InitializeLogSystem(InputFileName);
   
       // Input File Interpreter ####################################################################
       //string InputFileName = argv[1];
       InputFileData InputData;
       bool Check_InputFile_Flag = ReadCSV__Input_File_Interpreter(InputFileName, &InputData);
       if(Check_InputFile_Flag != true){
           LogSystem::write(" ->-> Error in input file!");
           LogSystem::write("");
           return 0;
       }
   
       // MagDataExplorer ###########################################################################
       MagDataProperties MagDataProp;
       bool Check_MagData_Flag;
       if(InputData.MagData_activate_flag){
           Check_MagData_Flag = MagData_Observer(InputData.MagDataPath, &MagDataProp);
           if(Check_MagData_Flag != true){
               LogSystem::write(" ->-> Error in MagData!");
               LogSystem::write("");
               LogSystem::write("");
               return 0;
           }
       }
   
       // NucDataExpolorer ###########################################################################
       NucDataProperties NucDataProp;
       bool Check_NucData_Flag;
       if(InputData.NucData_activate_flag){
           Check_NucData_Flag = NucData_Observer(InputData.NucDataPath, &NucDataProp);
           if(Check_NucData_Flag != true){
               LogSystem::write(" ->->Error in NucData!");
               LogSystem::write("");
               LogSystem::write("");
               return 0;
           }
       }
   
       // StructDataExplorer ########################################################################
       StructDataProperties StructDataProp;
       bool Check_StructData_Flag;
       if(InputData.StructData_activate_flag){
           Check_StructData_Flag = StructData_Observer(InputData.StructDataFilename, &StructDataProp);
           if(Check_StructData_Flag != true){
               LogSystem::write(" ->-> Error in StructData!");
               LogSystem::write("");
               LogSystem::write("");
           }
       }
   
       // Check Consistency of InputData and MagDataProp ############################################
   
   
   
   
   
       // Create Output Folder #####################################################################
       //mkdir((InputData.SANSDataFoldername + "/").c_str(), 0777);
       //LogSystem::initLog(InputData.SANSDataFoldername + "/NuMagSANSlog.txt");
       //LogSystem::write("Hello World!");
   
   
       // Start calculation based on loop modus or user selection ###################################
       //mkdir((InputData.SANSDataFoldername + "/").c_str(), 0777);    
   
       int Data_File_Index;
       size_t free_bytes, total_bytes;
       double used_mb1, used_mb2;
       double cummulated_mb = 0.0;
       if(InputData.Loop_Modus){
           LogSystem::write("loop modus active...");
           for(int k = InputData.Loop_From; k <= InputData.Loop_To; k++){
               Data_File_Index = k;
   
               cudaMemGetInfo(&free_bytes, &total_bytes);
               used_mb1 = (total_bytes - free_bytes) / 1024.0 / 1024.0;
   
               NuMagSANS_Calculator(&InputData, &NucDataProp, &MagDataProp, &StructDataProp, Data_File_Index);
   
               cudaMemGetInfo(&free_bytes, &total_bytes);
               used_mb2 = (total_bytes - free_bytes) / 1024.0 / 1024.0;
               cummulated_mb += used_mb2 - used_mb1;
   
               LogSystem::write("GPU Memory Check: cummulated bytes: " + std::to_string(cummulated_mb) + " MB, free bytes: " + std::to_string(free_bytes / 1024.0 /1024.0) + " MB");
   
           }
       }else{
           LogSystem::write("user selecting active...");
           for(int k = 0; k < InputData.Number_Of_User_Selections; k++){
               Data_File_Index = InputData.User_Selection_IndexArray[k];
   
               cudaMemGetInfo(&free_bytes, &total_bytes);
               used_mb1 = (total_bytes - free_bytes) / 1024.0 / 1024.0;
   
               NuMagSANS_Calculator(&InputData, &NucDataProp, &MagDataProp, &StructDataProp, Data_File_Index);
   
               cudaMemGetInfo(&free_bytes, &total_bytes);
               used_mb2 = (total_bytes - free_bytes) / 1024.0 / 1024.0;
               cummulated_mb += used_mb2 - used_mb1;
   
               LogSystem::write("GPU Memory Check: cummulated bytes: " + std::to_string(cummulated_mb) + " MB, free bytes: " + std::to_string(free_bytes / 1024.0 /1024.0) + " MB");
   
           }
       }
       
   
       LogSystem::close();
   
       return 0;
       
   }
   
   
