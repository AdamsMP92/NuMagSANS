// File         : MagSANSresponse.cu
// Author       : Michael Philipp ADAMS, M.Sc. 
// Company      : University of Luxembourg
// Department   : Department of Physics and Materials Sciences
// Group        : NanoMagnetism Group
// Group Leader : Prof. Andreas Michels
// Version      : 19 November 2024
// OS           : Linux Ubuntu
// Language     : CUDA C++

#include "NuMagSANSlib.h"

using namespace std;

int main(int argc, char* argv[]){
		
	// Input File Interpreter ####################################################################
	string InputFileName = argv[1];
	InputFileData InputData;
    bool Check_InputFile_Flag =	ReadCSV__Input_File_Interpreter(InputFileName, &InputData);
	if(Check_InputFile_Flag != true){
		std::cout << " ->-> Error in input file!" << "\n\n";
		return 0;
	}

	// MagDataObserver ###########################################################################
	MagDataProperties MagDataProp;
	bool Check_MagData_Flag;
	if(InputData.MagData_activate_flag){
		Check_MagData_Flag = MagData_Observer(InputData.MagDataPath, &MagDataProp);
		if(Check_MagData_Flag != true){
			std::cout << " ->-> Error in magnetization data!" << "\n\n";
			return 0;
		}
	}

	// NucDataObserver ###########################################################################
	NucDataProperties NucDataProp;
	bool Check_NucData_Flag;
	if(InputData.NucData_activate_flag){
		Check_NucData_Flag = NucData_Observer(InputData.NucDataPath, &NucDataProp);
		if(Check_NucData_Flag != true){
			std::cout << " ->-> Error in nuclear data!" << "\n\n";
			return 0;
		}
	}

	// StructDataObserver ########################################################################
	StructDataProperties StructDataProp;
	bool Check_StructData_Flag;
	if(InputData.StructData_activate_flag){
		Check_StructData_Flag = StructData_Observer(InputData.StructDataFilename, &StructDataProp);
		if(Check_StructData_Flag != true){
			std::cout << " ->-> Error in structure data!" << "\n\n";
		}
	}

	// Check Consistency of InputData and MagDataProp ############################################






	
	// Start calculation based on loop modus or user selection ###################################
	mkdir((InputData.SANSDataFoldername + "/").c_str(), 0777);	
	int Data_File_Index;
	if(InputData.Loop_Modus){
		std::cout << "Loop Modus active" << "\n";
		for(int k = InputData.Loop_From; k <= InputData.Loop_To; k++){
			Data_File_Index = k;
			NuMagSANS_Calculator(&InputData, &NucDataProp, &MagDataProp, &StructDataProp, Data_File_Index);
		}
	}else{
		std::cout << "User Selection active" << "\n";
		for(int k = 0; k < InputData.Number_Of_User_Selections; k++){
			Data_File_Index = InputData.User_Selection_IndexArray[k];
			NuMagSANS_Calculator(&InputData, &NucDataProp, &MagDataProp, &StructDataProp, Data_File_Index);
		}
	}
	
	return 0;
	
}


