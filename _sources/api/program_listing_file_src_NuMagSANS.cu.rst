
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
   
   int main(int argc, char* argv[]) {
   
       // Input File Name
       string InputFileName = argv[1];
   
       // Initialize LogFileSystem
       InitializeLogSystem(InputFileName);
   
       // Input File Interpreter ####################################################################
       InputFileData InputData;
       bool Check_InputFile_Flag = ReadCSV__Input_File_Interpreter(InputFileName, &InputData);
       if (Check_InputFile_Flag != true) {
           LogSystem::write(" ->-> Error in input file!");
           LogSystem::write("");
           return 0;
       }
   
       // MagDataExplorer ###########################################################################
       MagDataProperties MagDataProp;
       bool Check_MagData_Flag;
       if (InputData.MagData_activate_flag) {
           Check_MagData_Flag = MagData_Observer(InputData.MagDataPath, &MagDataProp, &InputData);
           if (Check_MagData_Flag != true) {
               LogSystem::write(" ->-> Error in MagData!");
               LogSystem::write("");
               LogSystem::write("");
               return 0;
           }
       }
   
       // NucDataExpolorer ###########################################################################
       NucDataProperties NucDataProp;
       bool Check_NucData_Flag;
       if (InputData.NucData_activate_flag) {
           Check_NucData_Flag = NucData_Observer(InputData.NucDataPath, &NucDataProp, &InputData);
           if (Check_NucData_Flag != true) {
               LogSystem::write(" ->->Error in NucData!");
               LogSystem::write("");
               LogSystem::write("");
               return 0;
           }
       }
   
       // StructDataExplorer ########################################################################
       StructDataProperties StructDataProp;
       bool Check_StructData_Flag;
       if (InputData.StructData_activate_flag) {
           if (InputData.StructDataLoop_flag) {
               Check_StructData_Flag =
                   StructDataLoop_Observer(InputData.StructDataPath, InputData.StructDataLoop_IndexArray, &StructDataProp);
           } else {
               Check_StructData_Flag = StructData_Observer(InputData.StructDataFilename, &StructDataProp);
           }
           if (Check_StructData_Flag != true) {
               LogSystem::write(" ->-> Error in StructData!");
               LogSystem::write("");
               LogSystem::write("");
               return 0;
           }
       }
   
       // RotDataExplorer ##########################################################################
       RotDataProperties RotDataProp;
       bool Check_RotData_Flag;
       if (InputData.RotData_activate_flag) {
           if (InputData.RotDataLoop_flag) {
               Check_RotData_Flag =
                   RotDataLoop_Observer(InputData.RotDataPath, InputData.RotDataLoop_IndexArray, &RotDataProp);
           } else {
               Check_RotData_Flag = RotData_Observer(InputData.RotDataFilename, &RotDataProp);
           }
           if (Check_RotData_Flag != true) {
               LogSystem::write(" ->-> Error in RotData!");
               LogSystem::write("");
               LogSystem::write("");
               return 0;
           }
       }
   
       if (InputData.StructData_activate_flag && InputData.RotData_activate_flag) {
           if (StructDataProp.Number_Of_Elements != RotDataProp.Number_Of_Elements) {
               LogSystem::write(" ->-> Error: StructData and RotData contain different numbers of entries!");
               LogSystem::write("StructData entries: " + std::to_string(StructDataProp.Number_Of_Elements));
               LogSystem::write("RotData entries: " + std::to_string(RotDataProp.Number_Of_Elements));
               LogSystem::write("");
               return 0;
           }
       }
   
       int EffectiveObjectCount = 0;
       if (InputData.MagData_activate_flag) {
           EffectiveObjectCount = MagDataProp.Number_Of_SubFolders;
       }
       if (InputData.NucData_activate_flag) {
           if (EffectiveObjectCount == 0) {
               EffectiveObjectCount = NucDataProp.Number_Of_SubFolders;
           } else if (EffectiveObjectCount != NucDataProp.Number_Of_SubFolders) {
               LogSystem::write(" ->-> Error: MagData and NucData contain different effective object counts!");
               LogSystem::write("MagData objects: " + std::to_string(MagDataProp.Number_Of_SubFolders));
               LogSystem::write("NucData objects: " + std::to_string(NucDataProp.Number_Of_SubFolders));
               LogSystem::write("");
               return 0;
           }
       }
       if (InputData.StructData_activate_flag && EffectiveObjectCount > 0) {
           if (StructDataProp.Number_Of_Elements != EffectiveObjectCount) {
               LogSystem::write(" ->-> Error: StructData entries do not match the effective object count!");
               LogSystem::write("StructData entries: " + std::to_string(StructDataProp.Number_Of_Elements));
               LogSystem::write("Effective object count: " + std::to_string(EffectiveObjectCount));
               LogSystem::write("");
               return 0;
           }
       }
       if (InputData.RotData_activate_flag && EffectiveObjectCount > 0) {
           if (RotDataProp.Number_Of_Elements != EffectiveObjectCount) {
               LogSystem::write(" ->-> Error: RotData entries do not match the effective object count!");
               LogSystem::write("RotData entries: " + std::to_string(RotDataProp.Number_Of_Elements));
               LogSystem::write("Effective object count: " + std::to_string(EffectiveObjectCount));
               LogSystem::write("");
               return 0;
           }
       }
   
       // Check Consistency of InputData and MagDataProp ############################################
   
       // Start calculation based on loop modus or user selection ###################################
       int Data_File_Index;
       if (InputData.Loop_Modus) {
           LogSystem::write("loop modus active...");
           for (int k = InputData.Loop_From; k <= InputData.Loop_To; k++) {
               Data_File_Index = k;
               NuMagSANS_Calculator(&InputData, &NucDataProp, &MagDataProp, &StructDataProp, &RotDataProp,
                                    Data_File_Index);
           }
       } else {
           LogSystem::write("user selecting active...");
           for (int k = 0; k < InputData.User_Selection_IndexArray.size(); k++) {
               Data_File_Index = InputData.User_Selection_IndexArray[k];
               NuMagSANS_Calculator(&InputData, &NucDataProp, &MagDataProp, &StructDataProp, &RotDataProp,
                                    Data_File_Index);
           }
       }
   
       LogSystem::close();
   
       return 0;
   }
