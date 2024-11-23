// File         : NuMagSANSlib_MagDataObserver.h
// Author       : Michael Philipp ADAMS, M.Sc.
// Company      : University of Luxembourg
// Department   : Department of Physics and Materials Sciences
// Group        : NanoMagnetism Group
// Group Leader : Prof. Andreas Michels
// Version      : 22 October 2024
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

struct MagDataProperties{

	string GlobalFolderPath;

	string SubFolderNames_Nom;
	int Number_Of_SubFolders;
	
	string SubFolder_FileNames_Nom;
    string SubFolder_FileNames_Type;
	int Number_Of_Files_In_SubFolder;

	int** NumberOfElements;
	int** NumberOfNonZeroMoments;
	
};

// ###############################################################################################################################################
// helper functions ##############################################################################################################################
// ###############################################################################################################################################

void get_GlobalMagDataPath(std::string Local_MagDataPath, MagDataProperties* MagDataProp){

    char tmp[PATH_MAX];
    getcwd(tmp, PATH_MAX);  // Get the current working directory
    std::string tmp_string = tmp;
	MagDataProp->GlobalFolderPath = tmp_string + "/" + Local_MagDataPath;
	cout << "found Global MagDataPath: " << MagDataProp->GlobalFolderPath << "\n\n";
    
}




bool check_number_of_elements_in_folders_MagData(MagDataProperties* MagDataProp){

	// Read number of elements in each subfolder
	int Number_Of_Elements[MagDataProp->Number_Of_SubFolders];
	string current_path;
	for(int k=0; k < MagDataProp->Number_Of_SubFolders; k++){
		current_path = MagDataProp->GlobalFolderPath + "/" + MagDataProp->SubFolderNames_Nom + "_" + std::to_string(k+1);
		Number_Of_Elements[k] = count_NumberOfElements(current_path);
	}
	cout << "\n";
			// check wether number of elements in each subfolder is the same
	bool Check_Flag = true;
	for(int k=1; k < MagDataProp->Number_Of_SubFolders; k++){
		if(Number_Of_Elements[k] != Number_Of_Elements[0]){
			Check_Flag = false;
			break;
		}
	}	
	if(Check_Flag){
		cout << "Each folder got the same number of elements: " << Number_Of_Elements[0] << "\n";
		MagDataProp->Number_Of_Files_In_SubFolder = Number_Of_Elements[0];
	}
	return Check_Flag;

}

bool check_element_names_MagData(MagDataProperties* MagDataProp){

	// read names of elements in folder with index 1
	std::string ElementNames[MagDataProp->Number_Of_Files_In_SubFolder];
	std::string	current_path;
	bool ElementNames_CheckFlag = false;
	std::string FileNames_Type[MagDataProp->Number_Of_SubFolders];
	std::string FileNames_Nom[MagDataProp->Number_Of_SubFolders];


	for(int k=0; k < MagDataProp->Number_Of_SubFolders; k++){

		ElementNames_CheckFlag = false;

		// read all element names in current_path directory
		current_path = MagDataProp->GlobalFolderPath + "/" + MagDataProp->SubFolderNames_Nom + "_" + std::to_string(k+1);
		read_FolderNames(current_path, ElementNames, MagDataProp->Number_Of_Files_In_SubFolder);

		// check element names in folder 1
		CheckStrings_ElementNames(ElementNames, MagDataProp->Number_Of_Files_In_SubFolder,\
									 &ElementNames_CheckFlag, &FileNames_Nom[k], &FileNames_Type[k]);

		if(ElementNames_CheckFlag != true){
			return false;
		}

	}

	// Compare Nom and Type	
	for(int k=1; k < MagDataProp->Number_Of_SubFolders; k++){
		if(FileNames_Type[k] != FileNames_Type[0] || FileNames_Nom[k] != FileNames_Nom[0]){
			return false;
		}
	}

	// Store correct Nom and Type
	MagDataProp->SubFolder_FileNames_Nom = FileNames_Nom[0];
	MagDataProp->SubFolder_FileNames_Type = FileNames_Type[0];

	return true;
	
}

bool check_Subfolders_MagData(MagDataProperties* MagDataProp){

	// count number of subfolders
	MagDataProp->Number_Of_SubFolders = count_NumberOfFolders(MagDataProp->GlobalFolderPath);
	cout << "Number of Subfolders: " << MagDataProp->Number_Of_SubFolders << "\n";
	
	// Read all subfolder names
	std::string SubFolderNames[MagDataProp->Number_Of_SubFolders];
	read_FolderNames(MagDataProp->GlobalFolderPath, SubFolderNames, MagDataProp->Number_Of_SubFolders);

	bool FolderNames_CheckFlag = false;
	CheckStrings_FolderNames(SubFolderNames, MagDataProp->Number_Of_SubFolders, &FolderNames_CheckFlag, &MagDataProp->SubFolderNames_Nom);

	cout << "We found the Nom: " << MagDataProp->SubFolderNames_Nom << "\n";
	cout << "FolderNames CheckFlag: " << FolderNames_CheckFlag << "\n\n";

	return FolderNames_CheckFlag;
	
}

bool check_Subfolder_FileNames_MagData(MagDataProperties* MagDataProp){

	bool Elements_CheckFlag = check_number_of_elements_in_folders_MagData(MagDataProp);
	
	if(Elements_CheckFlag){

		bool Element_Names_CheckFlag = check_element_names_MagData(MagDataProp);

		if(Element_Names_CheckFlag){
			cout << "We found the Nom: " << MagDataProp->SubFolder_FileNames_Nom << ", and the file Type: " << MagDataProp->SubFolder_FileNames_Type << "\n\n";
			return true;
		}
		else{
			return false;
		}
	}
	else{
		return false;
	}	
}

bool check_FileDimensions_MagData(MagDataProperties* MagDataProp){

	string filename;

	// allocate memory
	MagDataProp->NumberOfElements = new int*[MagDataProp->Number_Of_SubFolders];
	MagDataProp->NumberOfNonZeroMoments = new int*[MagDataProp->Number_Of_SubFolders];
	for(int i = 0; i < MagDataProp->Number_Of_SubFolders; i++){
		MagDataProp->NumberOfElements[i] = new int[MagDataProp->Number_Of_Files_In_SubFolder];
		MagDataProp->NumberOfNonZeroMoments[i] = new int[MagDataProp->Number_Of_Files_In_SubFolder];
	} 

	// read the file information
	for(int i = 0; i < MagDataProp->Number_Of_SubFolders; i++){
		for(int j = 0; j < MagDataProp->Number_Of_Files_In_SubFolder; j++){		
		
			filename = MagDataProp->GlobalFolderPath + "/" + MagDataProp->SubFolderNames_Nom + "_" + std::to_string(i+1) + "/" \
	    		     + MagDataProp->SubFolder_FileNames_Nom + "_" + std::to_string(j+1)  + "." + MagDataProp->SubFolder_FileNames_Type;

			NumberOfNonZeroMagneticMomentsInFile(&MagDataProp->NumberOfNonZeroMoments[i][j], &MagDataProp->NumberOfElements[i][j], filename);

		}
	}

	return true;
		
}

// Routine that checks number of subfolders in MagData directory
bool MagData_Observer(std::string Local_MagDataPath, MagDataProperties*MagDataProp){

	cout << "##########################################################################################" << "\n";
	cout << "## Run - MagData Directory Explorer ######################################################" << "\n";
	cout << "##########################################################################################" << "\n\n";

	bool CheckFlag = false;

	// get global path of the MagData directory
	get_GlobalMagDataPath(Local_MagDataPath, MagDataProp);

	// check the subfolder names
	bool Subfolder_CheckFlag = check_Subfolders_MagData(MagDataProp);

	// Check the files in each subfolder
	bool Subfolder_Elements_CheckFlag = check_Subfolder_FileNames_MagData(MagDataProp);


	bool FileDimensions_CheckFlag = check_FileDimensions_MagData(MagDataProp);


	cout << "Number Of Elements: " << MagDataProp->NumberOfElements[0][0] << "\n";
	cout << "Non-Zero magnetic moments: " << MagDataProp->NumberOfNonZeroMoments[0][0] << "\n";
	cout << "\n";

	cout << "##########################################################################################" << "\n";
	cout << "## Stop - MagData Directory Explorer #####################################################" << "\n";
	cout << "##########################################################################################" << "\n\n";

	if(Subfolder_CheckFlag && Subfolder_Elements_CheckFlag && FileDimensions_CheckFlag){
		CheckFlag = true;
	}

	return CheckFlag;

}

