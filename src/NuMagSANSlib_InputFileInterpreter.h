// File         : NuMagSANSlib_InputFileInterpreter.h
// Author       : Michael Philipp ADAMS, M.Sc.
// Company      : University of Luxembourg
// Department   : Department of Physics and Materials Sciences
// Group        : NanoMagnetism Group
// Group Leader : Prof. Andreas Michels
// Version      : 28 November 2024
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


// ##############################################################################################################################################
// ##############################################################################################################################################


struct InputFileData{

	string NucDataPath;
	string MagDataPath;
	string StructDataFilename;
	string SANSDataFoldername;
	string Fourier_Approach;

	bool Loop_Modus;
	int Loop_From;
	int Loop_To;
	
	string User_Selection;
	int Number_Of_User_Selections;
	int *User_Selection_IndexArray;

    int N_q; 
    int N_theta;
    int N_r; 
    int N_alpha;
    float q_max;
    float r_max;

    float Scattering_Volume_V;
	float cell_nuclear_sld;
    float cell_magnetization;
	float cuboid_cell_size_x;
	float cuboid_cell_size_y;
	float cuboid_cell_size_z;
    
    float RotMat_alpha;
    float RotMat_beta;
    float RotMat[9];
    float XYZ_Unit_Factor;

	float Polarization[3];

	bool NucData_activate_flag; 
	bool MagData_activate_flag;
	bool StructData_activate_flag;
	bool ExcludeZeroMoments_flag;
    
    bool Check_InputFile_Flag;
    
    bool output_fourier_correlation_matrix_flag;
    bool output_unpolarized_nuclear_SANS_cross_section_2D_flag;
    bool output_unpolarized_magnetic_SANS_cross_section_2D_flag;
	bool output_polarized_magnetic_SANS_cross_section_2D_flag;
	bool output_nuclear_magnetic_SANS_cross_section_2D_flag;
    bool output_spin_flip_magnetic_SANS_cross_section_2D_flag;
    bool output_chiral_magnetic_SANS_cross_section_2D_flag;
    bool output_pm_spin_flip_magnetic_SANS_cross_section_2D_flag;
    bool output_mp_spin_flip_magnetic_SANS_cross_section_2D_flag;
	bool output_pp_non_spin_flip_magnetic_SANS_cross_section_2D_flag;
	bool output_mm_non_spin_flip_magnetic_SANS_cross_section_2D_flag;
	bool output_p_sanspol_magnetic_SANS_cross_section_2D_flag;
	bool output_m_sanspol_magnetic_SANS_cross_section_2D_flag;

	bool output_unpolarized_nuclear_SANS_cross_section_1D_flag;
    bool output_unpolarized_magnetic_SANS_cross_section_1D_flag;
	bool output_polarized_magnetic_SANS_cross_section_1D_flag;
	bool output_nuclear_magnetic_SANS_cross_section_1D_flag;
    bool output_spin_flip_magnetic_SANS_cross_section_1D_flag;
    bool output_chiral_magnetic_SANS_cross_section_1D_flag;
    bool output_pm_spin_flip_SANS_cross_section_1D_flag;
    bool output_mp_spin_flip_SANS_cross_section_1D_flag;
    bool output_pp_non_spin_flip_SANS_cross_section_1D_flag;
    bool output_mm_non_spin_flip_SANS_cross_section_1D_flag;
    bool output_m_sanspol_cross_section_1D_flag; 
    bool output_p_sanspol_cross_section_1D_flag;

	bool output_nuclear_pair_distance_distribution_1D_flag;
	bool output_unpolarized_pair_distance_distribution_1D_flag;
	bool output_polarized_pair_distance_distribution_1D_flag;
	bool output_nuclear_magnetic_pair_distance_distribution_1D_flag;
	bool output_spin_flip_pair_distance_distribution_1D_flag;
	bool output_chiral_pair_distance_distribution_1D_flag;
	bool output_pm_spin_flip_pair_distance_distribution_1D_flag;
	bool output_mp_spin_flip_pair_distance_distribution_1D_flag;
	bool output_pp_non_spin_flip_pair_distance_distribution_1D_flag;
	bool output_mm_non_spin_flip_pair_distance_distribution_1D_flag;
	bool output_p_sanspol_pair_distance_distribution_1D_flag;
	bool output_m_sanspol_pair_distance_distribution_1D_flag;
	
	bool output_nuclear_correlation_function_1D_flag;
    bool output_unpolarized_correlation_function_1D_flag; 
    bool output_polarized_correlation_function_1D_flag;
    bool output_nuclear_magnetic_correlation_function_1D_flag;
    bool output_spin_flip_correlation_function_1D_flag;
    bool output_chiral_correlation_function_1D_flag;
    bool output_pm_spin_flip_correlation_function_1D_flag;
   	bool output_mp_spin_flip_correlation_function_1D_flag;
   	bool output_pp_non_spin_flip_correlation_function_1D_flag;
   	bool output_mm_non_spin_flip_correlation_function_1D_flag;
   	bool output_p_sanspol_correlation_function_1D_flag;
   	bool output_m_sanspol_correlation_function_1D_flag;

	bool output_nuclear_correlation_function_2D_flag;
    bool output_unpolarized_correlation_function_2D_flag;
    bool output_polarized_correlation_function_2D_flag; 
    bool output_nuclear_magnetic_correlation_function_2D_flag;
    bool output_spin_flip_correlation_function_2D_flag;
    bool output_chiral_correlation_function_2D_flag;
    bool output_pm_spin_flip_correlation_function_2D_flag; 
    bool output_mp_spin_flip_correlation_function_2D_flag;
    bool output_pp_non_spin_flip_correlation_function_2D_flag;
    bool output_mm_non_spin_flip_correlation_function_2D_flag;
    bool output_m_sanspol_correlation_function_2D_flag;
    bool output_p_sanspol_correlation_function_2D_flag;

};

// zy-rotation matrix where alpha is the polar angle and beta is the azimuth angle
void Compute_RotMat(float alpha, float beta, float* RotMat){

	// ordering of the rotation matrix:
	// RotMat = [RotMat[0], RotMat[3], RotMat[6]; ...
	//           RotMat[1], RotMat[4], RotMat[7]; ...
	//           RotMat[2], RotMat[5], RotMat[8]]


	 alpha = alpha * M_PI/180.0;
	 beta = beta * M_PI/180.0;

     RotMat[0] = cos(alpha) * cos(beta);
     RotMat[1] = cos(alpha) * sin(beta);
     RotMat[2] = -sin(alpha);
     RotMat[3] = -sin(beta);
     RotMat[4] = cos(beta);
     RotMat[5] = 0;
     RotMat[6] = cos(beta) * sin(alpha);
     RotMat[7] = sin(alpha) * sin(beta);
     RotMat[8] = cos(alpha);

 }


string extract_value(const string &line){
	size_t equality_symbol_index = line.find("=");
    size_t semicolon_symbol_index = line.find(";");
	size_t space_symbol_index = 0;
	string string_to_be_interpreted = line.substr(equality_symbol_index+1, semicolon_symbol_index-equality_symbol_index-1);
	while(string_to_be_interpreted.find(" ") <= string_to_be_interpreted.size()){
                                space_symbol_index = string_to_be_interpreted.find(" ");
                                string_to_be_interpreted.erase(space_symbol_index, 1);
        }
	return string_to_be_interpreted;
}

void parseBool(const std::string& line, const std::string& property, bool& variable, bool& flag){
	size_t found = line.find(property);
	if (found != std::string::npos) {
   		variable = (bool) stoi(extract_value(line));
		flag = true;
		cout << " -> " << property << " : " << variable << "\n"; 
    }
}

void parseInt(const std::string& line, const std::string& property, int& variable, bool& flag){
	size_t found = line.find(property);
	if (found != std::string::npos) {
   		variable = stoi(extract_value(line));
		flag = true;
		cout << " -> " << property << " : " << variable << "\n";
    }
}

void parseFloat(const std::string& line, const std::string& property, float& variable, bool& flag){
	size_t found = line.find(property);
	if (found != std::string::npos) {
   		variable = stof(extract_value(line));
		flag = true;
		cout << " -> " << property << " : " << variable << "\n";
    }
}

void parseString(const std::string& line, const std::string& property, std::string& variable, bool& flag){
	size_t found = line.find(property);
	if (found != std::string::npos) {
   		variable = extract_value(line);
		flag = true;
		cout << " -> " << property << " : " << variable << "\n";
    }
}

 // Interpreting User_Selection to integer array ##########################################
void User_Selection_To_Int_Array(int *Idx, int Number_of_Comma, string my_string){
 
     int Symbol_Idx[Number_of_Comma + 2];
     string Symbol_Start = "{";
     string Symbol_End = "}";
     string Symbol_Comma = ",";
     int Symbol_Counter = 0;
     for(int i = 0; i < my_string.size(); i++){
         if(my_string.at(i) == Symbol_Start.at(0) || my_string.at(i) == Symbol_End.at(0) || my_string.at(i) == Symbol_Comma.at(0)){
             Symbol_Idx[Symbol_Counter] = i;
             Symbol_Counter += 1;
         }
     }

     string Sub_String;
     int space_symbol_index;
     for(int i = 0; i <Number_of_Comma + 1; i++){
         Sub_String = my_string.substr(Symbol_Idx[i]+1, Symbol_Idx[i+1]-Symbol_Idx[i]-1);
         while(Sub_String.find(" ") <= Sub_String.size()){
             space_symbol_index = Sub_String.find(" ");
             Sub_String.erase(space_symbol_index, 1);
         }
         Idx[i] = stoi(Sub_String);
     }
}

// Find Number of Comma in string #########################################################
void Number_Of_Comma_In_String(int *N, string my_string){
     int String_Length = my_string.size();
     int Comma_Counter = 0;
     string Comma = ",";

     for(int i = 0; i < String_Length; i++){
         if(my_string.at(i) == Comma.at(0)){
             Comma_Counter += 1;
         }
     }
     *N = Comma_Counter;
}


bool ReadCSV__Input_File_Interpreter(string filename, InputFileData*InputData){

	cout << "##########################################################################################" << "\n";
	cout << "## Run - Input File Interpreter ##########################################################" << "\n";
	cout << "##########################################################################################" << "\n\n";

	bool Check_InputFile_Flag = false;

	ifstream fin;
	fin.open(filename);
	string line;
	string string_to_be_interpreted;

	int line_counter = 0;

	bool Check_Flag[92];
	for(int k=0; k<92; k++){
		Check_Flag[k] = false;
	}

	cout << "Interpreter found the following commands:" << "\n\n";

	while(getline(fin, line)){
		line_counter += 1;

		parseBool(line, "Loop_Modus", InputData->Loop_Modus, Check_Flag[0]);
    	parseInt(line, "Loop_From", InputData->Loop_From, Check_Flag[1]);
		parseInt(line, "Loop_To", InputData->Loop_To, Check_Flag[2]);		
		parseInt(line, "Number_Of_q_Points", InputData->N_q, Check_Flag[3]);
		parseInt(line, "Number_Of_theta_Points", InputData->N_theta, Check_Flag[4]);
		parseInt(line, "Number_Of_r_Points", InputData->N_r, Check_Flag[5]);
		parseInt(line, "Number_Of_alpha_Points", InputData->N_alpha, Check_Flag[6]);
		parseString(line, "User_Selection", InputData->User_Selection, Check_Flag[7]);
		parseFloat(line, "q_max", InputData->q_max, Check_Flag[8]);
		parseFloat(line, "r_max", InputData->r_max, Check_Flag[9]);
		parseFloat(line, "XYZ_Unit_Fact", InputData->XYZ_Unit_Factor, Check_Flag[10]);
		parseFloat(line, "RotMat_alpha", InputData->RotMat_alpha, Check_Flag[11]);
		parseFloat(line, "RotMat_beta", InputData->RotMat_beta, Check_Flag[12]);
		parseBool(line, "Fourier_Gamma", InputData->output_fourier_correlation_matrix_flag, Check_Flag[13]);
		parseBool(line, "Unpolarized_2D", InputData->output_unpolarized_magnetic_SANS_cross_section_2D_flag, Check_Flag[14]);
		parseBool(line, "AVG_SpinFlip_2D", InputData->output_spin_flip_magnetic_SANS_cross_section_2D_flag, Check_Flag[15]);
		parseBool(line, "Chiral_2D", InputData->output_chiral_magnetic_SANS_cross_section_2D_flag, Check_Flag[16]);
		parseBool(line, "PM_SpinFlip_2D", InputData->output_pm_spin_flip_magnetic_SANS_cross_section_2D_flag, Check_Flag[17]);
		parseBool(line, "MP_SpinFlip_2D", InputData->output_mp_spin_flip_magnetic_SANS_cross_section_2D_flag, Check_Flag[18]);
		parseBool(line, "Unpolarized_1D", InputData->output_unpolarized_magnetic_SANS_cross_section_1D_flag, Check_Flag[19]);
		parseBool(line, "AVG_SpinFlip_1D", InputData->output_spin_flip_magnetic_SANS_cross_section_1D_flag, Check_Flag[20]);
		parseBool(line, "Chiral_1D", InputData->output_chiral_magnetic_SANS_cross_section_1D_flag, Check_Flag[21]);
		parseBool(line, "Unpolarized_PairDist_1D", InputData->output_unpolarized_pair_distance_distribution_1D_flag, Check_Flag[22]);
		parseBool(line, "Unpolarized_Corr_1D", InputData->output_unpolarized_correlation_function_1D_flag, Check_Flag[23]);
		parseBool(line, "SpinFlip_PairDist_1D", InputData->output_spin_flip_pair_distance_distribution_1D_flag, Check_Flag[24]);
		parseBool(line, "SpinFlip_Corr_1D", InputData->output_spin_flip_correlation_function_1D_flag, Check_Flag[25]);
		parseBool(line, "Chiral_PairDist_1D", InputData->output_chiral_pair_distance_distribution_1D_flag, Check_Flag[26]);
		parseBool(line, "Chiral_Corr_1D", InputData->output_chiral_correlation_function_1D_flag, Check_Flag[27]);
		parseBool(line, "Unpolarized_Corr_2D", InputData->output_unpolarized_correlation_function_2D_flag, Check_Flag[28]);
		parseBool(line, "SpinFlip_Corr_2D", InputData->output_spin_flip_correlation_function_2D_flag, Check_Flag[29]);
		parseBool(line, "Chiral_Corr_2D", InputData->output_chiral_correlation_function_2D_flag, Check_Flag[30]);
		parseString(line, "MagDataPath", InputData->MagDataPath, Check_Flag[31]);
		parseString(line, "foldernameSANSData", InputData->SANSDataFoldername, Check_Flag[32]);
		parseFloat(line, "Scattering_Volume_V", InputData->Scattering_Volume_V, Check_Flag[33]);
		parseString(line, "Fourier_Approach", InputData->Fourier_Approach, Check_Flag[34]);
		parseFloat(line, "Cell_Magnetization", InputData->cell_magnetization, Check_Flag[35]);
		parseFloat(line, "Cuboid_Cell_Size_x", InputData->cuboid_cell_size_x, Check_Flag[36]);
		parseFloat(line, "Cuboid_Cell_Size_y", InputData->cuboid_cell_size_y, Check_Flag[37]);
		parseFloat(line, "Cuboid_Cell_Size_z", InputData->cuboid_cell_size_z, Check_Flag[38]);
		parseFloat(line, "Polarization_x", InputData->Polarization[0], Check_Flag[39]);
		parseFloat(line, "Polarization_y", InputData->Polarization[1], Check_Flag[40]);
		parseFloat(line, "Polarization_z", InputData->Polarization[2], Check_Flag[41]);
		parseBool(line, "Polarized_2D", InputData->output_polarized_magnetic_SANS_cross_section_2D_flag, Check_Flag[42]);
		parseBool(line, "NuclearMagnetic_2D", InputData->output_nuclear_magnetic_SANS_cross_section_2D_flag, Check_Flag[43]);
		parseBool(line, "PP_NonSpinFlip_2D", InputData->output_pp_non_spin_flip_magnetic_SANS_cross_section_2D_flag, Check_Flag[44]);
		parseBool(line, "MM_NonSpinFlip_2D", InputData->output_mm_non_spin_flip_magnetic_SANS_cross_section_2D_flag, Check_Flag[45]);
		parseBool(line, "P_SANSPOL_2D", InputData->output_p_sanspol_magnetic_SANS_cross_section_2D_flag, Check_Flag[46]);
		parseBool(line, "M_SANSPOL_2D", InputData->output_m_sanspol_magnetic_SANS_cross_section_2D_flag, Check_Flag[47]);
		parseString(line, "NucDataPath", InputData->NucDataPath, Check_Flag[48]);
		parseBool(line, "NucData_activate", InputData->NucData_activate_flag, Check_Flag[49]);
		parseBool(line, "MagData_activate", InputData->MagData_activate_flag, Check_Flag[50]);
		parseBool(line, "Exclude_Zero_Moments", InputData->ExcludeZeroMoments_flag, Check_Flag[51]);
		parseBool(line, "Nuclear_1D", InputData->output_unpolarized_nuclear_SANS_cross_section_1D_flag, Check_Flag[52]);
		parseString(line, "StructDataFilename", InputData->StructDataFilename, Check_Flag[53]);
		parseBool(line, "StructData_activate", InputData->StructData_activate_flag, Check_Flag[54]);
		parseFloat(line, "Cell_Nuclear_SLD", InputData->cell_nuclear_sld, Check_Flag[55]);
		parseBool(line, "Polarized_1D", InputData->output_polarized_magnetic_SANS_cross_section_1D_flag, Check_Flag[56]);
		parseBool(line, "NuclearMagnetic_1D", InputData->output_nuclear_magnetic_SANS_cross_section_1D_flag, Check_Flag[57]);
		parseBool(line, "Polarized_PairDist_1D", InputData->output_polarized_pair_distance_distribution_1D_flag, Check_Flag[58]);
		parseBool(line, "Polarized_Corr_1D", InputData->output_polarized_correlation_function_1D_flag, Check_Flag[59]);
		parseBool(line, "Nuclear_PairDist_1D", InputData->output_nuclear_pair_distance_distribution_1D_flag, Check_Flag[60]);
		parseBool(line, "Nuclear_Corr_1D", InputData->output_nuclear_correlation_function_1D_flag, Check_Flag[61]);
		parseBool(line, "PM_SpinFlip_1D", InputData->output_pm_spin_flip_SANS_cross_section_1D_flag, Check_Flag[62]);
		parseBool(line, "MP_SpinFlip_1D", InputData->output_mp_spin_flip_SANS_cross_section_1D_flag, Check_Flag[63]);
		parseBool(line, "PP_NonSpinFlip_1D", InputData->output_pp_non_spin_flip_SANS_cross_section_1D_flag, Check_Flag[64]);
		parseBool(line, "MM_NonSpinFlip_1D", InputData->output_mm_non_spin_flip_SANS_cross_section_1D_flag, Check_Flag[65]);
		parseBool(line, "P_SANSPOL", InputData->output_p_sanspol_cross_section_1D_flag, Check_Flag[66]);
		parseBool(line, "M_SANSPOL", InputData->output_m_sanspol_cross_section_1D_flag, Check_Flag[67]);
		parseBool(line, "Nuclear_Corr_2D", InputData->output_nuclear_correlation_function_2D_flag, Check_Flag[68]);
		parseBool(line, "Polarized_Corr_2D", InputData->output_polarized_correlation_function_2D_flag, Check_Flag[69]);
		parseBool(line, "NuclearMagnetic_Corr_2D", InputData->output_nuclear_magnetic_correlation_function_2D_flag, Check_Flag[70]);
		parseBool(line, "PM_SpinFlip_Corr_2D", InputData->output_pm_spin_flip_correlation_function_2D_flag, Check_Flag[71]);
		parseBool(line, "MP_SpinFlip_Corr_2D", InputData->output_pm_spin_flip_correlation_function_2D_flag, Check_Flag[72]);
		parseBool(line, "PP_NonSpinFlip_Corr_2D", InputData->output_pp_non_spin_flip_correlation_function_2D_flag, Check_Flag[73]);
		parseBool(line, "MM_NonSpinFlip_Corr_2D", InputData->output_mm_non_spin_flip_correlation_function_2D_flag, Check_Flag[74]);
		parseBool(line, "P_SANSPOL_2D", InputData->output_p_sanspol_correlation_function_2D_flag, Check_Flag[75]);
		parseBool(line, "M_SANSPOL_2D", InputData->output_m_sanspol_correlation_function_2D_flag, Check_Flag[76]);
		parseBool(line, "NucMag_PairDist_1D", InputData->output_nuclear_magnetic_pair_distance_distribution_1D_flag, Check_Flag[77]);
		parseBool(line, "PM_SpinFlip_PairDist_1D", InputData->output_pm_spin_flip_pair_distance_distribution_1D_flag, Check_Flag[78]);
		parseBool(line, "MP_SpinFlip_PairDist_1D", InputData->output_mp_spin_flip_pair_distance_distribution_1D_flag, Check_Flag[79]);
		parseBool(line, "PP_NonSpinFlip_PairDist_1D", InputData->output_pp_non_spin_flip_pair_distance_distribution_1D_flag, Check_Flag[80]);
		parseBool(line, "MM_NonSpinFlip_PairDist_1D", InputData->output_mm_non_spin_flip_pair_distance_distribution_1D_flag, Check_Flag[81]);
		parseBool(line, "P_SANSPOL_PairDist_1D", InputData->output_p_sanspol_pair_distance_distribution_1D_flag, Check_Flag[82]);
		parseBool(line, "M_SANSPOL_PairDist_1D", InputData->output_m_sanspol_pair_distance_distribution_1D_flag, Check_Flag[83]);
		parseBool(line, "NucMag_Corr_1D", InputData->output_nuclear_magnetic_correlation_function_1D_flag, Check_Flag[84]);
		parseBool(line, "PM_SpinFlip_Corr_1D", InputData->output_pm_spin_flip_correlation_function_1D_flag, Check_Flag[85]);
		parseBool(line, "MP_SpinFlip_Corr_1D", InputData->output_mp_spin_flip_correlation_function_1D_flag, Check_Flag[86]);
		parseBool(line, "PP_NonSpinFlip_Corr_1D", InputData->output_pp_non_spin_flip_correlation_function_1D_flag, Check_Flag[87]);
		parseBool(line, "MM_NonSpinFlip_Corr_1D", InputData->output_mm_non_spin_flip_correlation_function_1D_flag, Check_Flag[88]);
		parseBool(line, "P_SANSPOL_Corr_1D", InputData->output_p_sanspol_correlation_function_1D_flag, Check_Flag[89]);
		parseBool(line, "M_SANSPOL_Corr_1D", InputData->output_m_sanspol_correlation_function_1D_flag, Check_Flag[90]);
		parseBool(line, "Nuclear_2D", InputData->output_unpolarized_nuclear_SANS_cross_section_2D_flag, Check_Flag[91]);


	}		


	fin.close();


	// Transfer the User_Selection to integer array
	if(Check_Flag[7]){
		int Number_of_Comma = 0;
		Number_Of_Comma_In_String(&Number_of_Comma, InputData->User_Selection);
		InputData->User_Selection_IndexArray = new int[Number_of_Comma+1];
		User_Selection_To_Int_Array(InputData->User_Selection_IndexArray, Number_of_Comma, InputData->User_Selection);

		cout << "\n\n";
		cout << "Check UserSelection entries that are transfered to integer array:" << "\n";
		for(int k = 0; k < Number_of_Comma+1; k++){
			cout << "  UserSelection: " << k << " : " << InputData->User_Selection_IndexArray[k] << "\n";
		}
		InputData->Number_Of_User_Selections = Number_of_Comma + 1;
	}


	// compute rotation matrix
	if(Check_Flag[11] && Check_Flag[12]){
		Compute_RotMat(InputData->RotMat_alpha, InputData->RotMat_beta, InputData->RotMat);
	}
	cout << "Rotation Matrix: " << "\n";
	cout << InputData->RotMat[0] << " " << InputData->RotMat[3] << " " << InputData->RotMat[6] << "\n";
	cout << InputData->RotMat[1] << " " << InputData->RotMat[4] << " " << InputData->RotMat[7] << "\n";
	cout << InputData->RotMat[2] << " " << InputData->RotMat[5] << " " << InputData->RotMat[8] << "\n";
	cout << "\n";
	

	cout << "\n\n";

	// Check the Error Flags
	bool Error_Detect = true;
	for(int k = 0; k < 92; k++){
		if(Check_Flag[k] != 1){
			cout << "Error Check Flag " << k << "\n";
			Error_Detect = false;
		}
	}

	// Check Polarization
	float P_norm = sqrt(pow(InputData->Polarization[0], 2) \
					  + pow(InputData->Polarization[1], 2) \
					  + pow(InputData->Polarization[2], 2));
	if(P_norm == 0){
		cout << "Error: Polarization magnitude is equal to zero!!" << "\n";
	}
	if(P_norm != 0 && P_norm != 1){
		InputData->Polarization[0] = (InputData->Polarization[0])/P_norm;
		InputData->Polarization[1] = (InputData->Polarization[1])/P_norm;
		InputData->Polarization[2] = (InputData->Polarization[2])/P_norm;
		cout << "Polarization normalized to 1!" << "\n";
	}



	// Give information on errors
	if(Error_Detect){
		cout << " ->-> No Errors Detected" << "\n";
		Check_InputFile_Flag = true;
	}

	cout << "\n\n";
	cout << "##########################################################################################" << "\n";
	cout << "## Stop - Input File Interpreter #########################################################" << "\n";
	cout << "##########################################################################################" << "\n\n";

	return Check_InputFile_Flag;

}
