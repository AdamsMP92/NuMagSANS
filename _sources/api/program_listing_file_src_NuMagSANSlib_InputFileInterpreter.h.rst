
.. _program_listing_file_src_NuMagSANSlib_InputFileInterpreter.h:

Program Listing for File NuMagSANSlib_InputFileInterpreter.h
============================================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_NuMagSANSlib_InputFileInterpreter.h>` (``src/NuMagSANSlib_InputFileInterpreter.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   // File         : NuMagSANSlib_InputFileInterpreter.h
   // Author       : Dr. Michael Philipp ADAMS
   // Company      : University of Luxembourg
   // Department   : Department of Physics and Materials Sciences
   // Group        : NanoMagnetism Group
   // Group Leader : Prof. Andreas Michels
   // Version      : 28 October 2025
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
   #include <algorithm>
   #include <cctype>
   
   using namespace std;
   
   
   // ##############################################################################################################################################
   // ##############################################################################################################################################
   
   
   struct OutFlag {
       bool Nuclear = false;
       bool Unpolarized = false;
       bool Polarized = false;
       bool NuclearMagnetic= false;
       bool SpinFlip= false;
       bool Chiral= false;
       bool PM_SpinFlip = false;
       bool MP_SpinFlip = false;
       bool PP_NonSpinFlip = false; 
       bool MM_NonSpinFlip = false;
       bool P_SANSPOL = false;
       bool M_SANSPOL = false;
   };
   
   struct OutFlagQuant {
       OutFlag SANS2D{};
       OutFlag SANS1D{};
       OutFlag Corr2D{};
       OutFlag Corr1D{};
       OutFlag PairDist1D{};
   };
   
   struct PrefixEntry {
       const char* name;
       bool OutFlag::* member;
   };
   
   struct SuffixEntry {
       const char* name;
       OutFlag OutFlagQuant::* quantity;
   };
   
   
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
       std::vector<int> User_Selection_IndexArray;
   
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
   
       int k_max; 
       bool AngularSpec_activate_flag;
       
       bool output_fourier_correlation_matrix_flag;
   
       OutFlagQuant OutFlags{};
   
   };
   
   
   inline bool any_active(const OutFlag& f)
   {
       return  f.Nuclear
            || f.Unpolarized
            || f.Polarized
            || f.NuclearMagnetic
            || f.SpinFlip
            || f.Chiral
            || f.PM_SpinFlip
            || f.MP_SpinFlip
            || f.PP_NonSpinFlip
            || f.MM_NonSpinFlip
            || f.P_SANSPOL
            || f.M_SANSPOL;
   }
   
   
   
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
   
    
   
   
   
   static std::string trim(const std::string& s)
   {
       size_t start = s.find_first_not_of(" \t\r\n");
       size_t end   = s.find_last_not_of(" \t\r\n");
   
       if (start == std::string::npos)
           return "";
   
       return s.substr(start, end - start + 1);
   }
   
   
   static bool try_extract_value(const std::string& line,
                                 const std::string& property,
                                 std::string& value)
   {
       // Remove leading and trailing whitespace from the entire line
       std::string trimmed = trim(line);
   
       // If the line is shorter than the property name,
       // it cannot possibly match
       if (trimmed.size() < property.size())
           return false;
   
       // The line must start with the exact property name.
       // This ensures prefix matching but does NOT yet validate
       // token boundaries.
       if (trimmed.compare(0, property.size(), property) != 0)
           return false;
   
       // Position right after the property name
       size_t pos = property.size();
   
       // Skip any whitespace between the property and the '=' symbol
       while (pos < trimmed.size() &&
              std::isspace(static_cast<unsigned char>(trimmed[pos])))
           ++pos;
   
       // The next non-whitespace character MUST be '='.
       // This enforces an exact key match and prevents
       // matches like "q_max_extra".
       if (pos >= trimmed.size() || trimmed[pos] != '=')
           return false;
   
       // Find the terminating semicolon after '='
       size_t semi_pos = trimmed.find(';', pos);
       if (semi_pos == std::string::npos)
           return false;
   
       // Extract the substring between '=' and ';'
       // and trim surrounding whitespace
       value = trim(trimmed.substr(pos + 1,
                                   semi_pos - pos - 1));
   
       return true;
   }
    
   
   void parseBool(const std::string& line,
                  const std::string& property,
                  bool& variable,
                  bool& flag)
   {
       std::string value;
   
       if (!try_extract_value(line, property, value))
           return;
   
       variable = static_cast<bool>(std::stoi(value));
       flag = true;
   
       LogSystem::write(" -> " + property +
                        " : " +
                        std::string(variable ? "true" : "false"));
   }
   
   void parseInt(const std::string& line,
                 const std::string& property,
                 int& variable,
                 bool& flag)
   {
       std::string value;
   
       if (!try_extract_value(line, property, value))
           return;
   
       variable = std::stoi(value);
       flag = true;
   
       LogSystem::write(" -> " + property +
                        " : " +
                        std::to_string(variable));
   }
   
   void parseFloat(const std::string& line,
                   const std::string& property,
                   float& variable,
                   bool& flag)
   {
       std::string value;
   
       if (!try_extract_value(line, property, value))
           return;
   
       variable = std::stof(value);
       flag = true;
   
       LogSystem::write(" -> " + property +
                        " : " +
                        std::to_string(variable));
   }
   
   void parseString(const std::string& line,
                    const std::string& property,
                    std::string& variable,
                    bool& flag)
   {
       if (!try_extract_value(line, property, variable))
           return;
   
       flag = true;
   
       LogSystem::write(" -> " + property +
                        " : " +
                        variable);
   }
   
   void parseStringNoCout(const std::string& line,
                          const std::string& property,
                          std::string& variable,
                          bool& flag)
   {
       if (!try_extract_value(line, property, variable))
           return;
   
       flag = true;
   }
   
   
   
   
   static std::vector<int> parse_int_list(const std::string& s)
   {
       std::vector<int> result;
   
       std::string trimmed = trim(s);
   
       if (trimmed.front() != '{' || trimmed.back() != '}')
           throw std::runtime_error("Invalid User_Selection format");
   
       std::string content = trimmed.substr(1, trimmed.size() - 2);
   
       std::stringstream ss(content);
       std::string item;
   
       while (std::getline(ss, item, ',')) {
           result.push_back(std::stoi(trim(item)));
       }
   
       return result;
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
   
   // Pre-Parsing for the initialization of the LogFile System
   void InitializeLogSystem(string filename){
   
       ifstream fin;
       fin.open(filename);
       string line;
       string string_to_be_interpreted;
       bool checkFlag = false;
   
       string SANSDataFoldername;
   
       while(getline(fin, line)){
           parseStringNoCout(line, "foldernameSANSData", SANSDataFoldername, checkFlag);
       }
       fin.close();
   
       if(checkFlag){
           mkdir((SANSDataFoldername + "/").c_str(), 0777);
               LogSystem::initLog(SANSDataFoldername + "/NuMagSANSlog.txt");
               //LogSystem::write("LogFile Initialized");
       }
       else {
               std::cerr << "Error: 'foldernameSANSData' not found in input file '"<< filename << "'. Logging not initialized.\n";
           std::cerr << "Aborting program.\n";
               std::exit(EXIT_FAILURE);
       }
   
   }
   
   
   
   
   
   
   // template for options that are parsed by the InputFileInterpreter
   template<typename T>
   struct Option {
       std::string key;
       T* target;
       bool required;
       bool found = false;
   };
   
   
   
   
   
   
   bool ReadCSV__Input_File_Interpreter(string filename, InputFileData*InputData){
   
   
       LogSystem::write("##########################################################################################");
       LogSystem::write("## Run - Input File Interpreter ##########################################################");
       LogSystem::write("##########################################################################################");
   
   
       ifstream fin(filename);
       if (!fin.is_open()) {
           LogSystem::write("Could not open input file.");
           return false;
       }
   
       // generalisation of the output options 
       std::vector<PrefixEntry> prefixes = {
           {"Nuclear",        &OutFlag::Nuclear},
           {"Unpolarized",    &OutFlag::Unpolarized},
           {"Polarized",      &OutFlag::Polarized},
           {"NuclearMagnetic",&OutFlag::NuclearMagnetic},
           {"SpinFlip",       &OutFlag::SpinFlip},
           {"Chiral",         &OutFlag::Chiral},
           {"PM_SpinFlip",    &OutFlag::PM_SpinFlip},
           {"MP_SpinFlip",    &OutFlag::MP_SpinFlip},
           {"PP_NonSpinFlip", &OutFlag::PP_NonSpinFlip},
           {"MM_NonSpinFlip", &OutFlag::MM_NonSpinFlip},
           {"P_SANSPOL",      &OutFlag::P_SANSPOL},
           {"M_SANSPOL",      &OutFlag::M_SANSPOL}
       };
   
       std::vector<SuffixEntry> suffixes = {
           {"_2D",            &OutFlagQuant::SANS2D},
           {"_1D",            &OutFlagQuant::SANS1D},
           {"_Corr_2D",       &OutFlagQuant::Corr2D},
           {"_Corr_1D",       &OutFlagQuant::Corr1D},
           {"_PairDist_1D",   &OutFlagQuant::PairDist1D}
       };
   
       std::vector<Option<bool>> Output_options;
       for (auto& p : prefixes)
       {
           for (auto& s : suffixes)
           {
               std::string key = std::string(p.name) + s.name;
               Output_options.push_back({
                   key,
                   &( (InputData->OutFlags.*(s.quantity)).*(p.member) ),
                   true
               }); 
           }
       }
   
       // bool options //////////////////////////////
       std::vector<Option<bool>> bool_options = {
   
           {"Loop_Modus", &InputData->Loop_Modus, true},
           {"NucData_activate", &InputData->NucData_activate_flag, true},
           {"MagData_activate", &InputData->MagData_activate_flag, true},
           {"StructData_activate", &InputData->StructData_activate_flag, true},
           {"Exclude_Zero_Moments", &InputData->ExcludeZeroMoments_flag, true},
           {"Angular_Spec", &InputData->AngularSpec_activate_flag, true},
           {"Fourier_Gamma", &InputData->output_fourier_correlation_matrix_flag, true}
   
       };
   
       // int options /////////////////////////////////
       std::vector<Option<int>> int_options = {
   
           {"Loop_From", &InputData->Loop_From, true},
           {"Loop_To", &InputData->Loop_To, true},
           {"Number_Of_q_Points", &InputData->N_q, true},
           {"Number_Of_theta_Points", &InputData->N_theta, true},
           {"Number_Of_r_Points", &InputData->N_r, true},
           {"Number_Of_alpha_Points", &InputData->N_alpha, true},
           {"k_max", &InputData->k_max, true}
   
       };
   
       // float options ///////////////////////////////
       std::vector<Option<float>> float_options = {
   
           {"q_max", &InputData->q_max, true},
           {"r_max", &InputData->r_max, true},
           {"Scattering_Volume_V", &InputData->Scattering_Volume_V, true},
           {"Cell_Nuclear_SLD", &InputData->cell_nuclear_sld, true},
           {"Cell_Magnetization", &InputData->cell_magnetization, true},
           {"Cuboid_Cell_Size_x", &InputData->cuboid_cell_size_x, true},
           {"Cuboid_Cell_Size_y", &InputData->cuboid_cell_size_y, true},
           {"Cuboid_Cell_Size_z", &InputData->cuboid_cell_size_z, true},
           {"RotMat_alpha", &InputData->RotMat_alpha, true},
           {"RotMat_beta", &InputData->RotMat_beta, true},
           {"XYZ_Unit_Factor", &InputData->XYZ_Unit_Factor, true},
           {"Polarization_x", &InputData->Polarization[0], true},
           {"Polarization_y", &InputData->Polarization[1], true},
           {"Polarization_z", &InputData->Polarization[2], true}
   
       };
   
       // string options
       std::vector<Option<std::string>> string_options = {
   
           {"NucDataPath", &InputData->NucDataPath, true},
           {"MagDataPath", &InputData->MagDataPath, true},
           {"StructDataFilename", &InputData->StructDataFilename, true},
           {"foldernameSANSData", &InputData->SANSDataFoldername, true},
           {"Fourier_Approach", &InputData->Fourier_Approach, true},
           {"User_Selection", &InputData->User_Selection, true}
   
       };
   
   
       std::string line;
       while(getline(fin, line)) {
           std::string value;
   
           for (auto& opt : Output_options) {
               if (try_extract_value(line, opt.key, value)) {
                   if (opt.found) {
                       LogSystem::write("Duplicate definition of option: " + opt.key);
                       return false;
                   }
                   *opt.target = static_cast<bool>(std::stoi(value));
                   opt.found = true;
                   LogSystem::write(" -> " + opt.key + " : " + (*opt.target ? "true" : "false"));
               }
           }
   
           for (auto& opt : bool_options) {
               if (try_extract_value(line, opt.key, value)) {
                   if (opt.found) {
                       LogSystem::write("Duplicate definition of option: " + opt.key);
                       return false;
                   }
                   *opt.target = static_cast<bool>(std::stoi(value));
                   opt.found = true;
                   LogSystem::write(" -> " + opt.key + " : " + (*opt.target ? "true" : "false"));
               }
           }
   
           for (auto& opt : int_options) {
               if (try_extract_value(line, opt.key, value)) {
                   if (opt.found) {
                       LogSystem::write("Duplicate definition of option: " + opt.key);
                       return false;
                   }
                   *opt.target = std::stoi(value);
                   opt.found = true;
                   LogSystem::write(" -> " + opt.key + " : " + std::to_string(*opt.target));
               }
           }
   
           for (auto& opt : float_options) {
               if (try_extract_value(line, opt.key, value)) {
                   if (opt.found) {
                       LogSystem::write("Duplicate definition of option: " + opt.key);
                       return false;
                   }
                   *opt.target = std::stof(value);
                   opt.found = true;
                   LogSystem::write(" -> " + opt.key + " : " + std::to_string(*opt.target));
               }
           }
   
           for (auto& opt : string_options) {
               if (try_extract_value(line, opt.key, value)) {
                   if (opt.found) {
                       LogSystem::write("Duplicate definition of option: " + opt.key);
                       return false;
                   }
                   *opt.target = value;
                   opt.found = true;
                   LogSystem::write(" -> " + opt.key + " : " + *opt.target);
               }
           }
   
       }
       fin.close();
   
   
   
   
       for (auto& opt : string_options)
       if (opt.key == "User_Selection" && opt.found)
       {
           InputData->User_Selection_IndexArray =
               parse_int_list(InputData->User_Selection);
   
           LogSystem::write("Check UserSelection entries that are transferred to integer array:");
   
           for (size_t k = 0; k < InputData->User_Selection_IndexArray.size(); ++k)
           {
               LogSystem::write(
                   " UserSelection: " +
                   std::to_string(k) + " : " +
                   std::to_string(InputData->User_Selection_IndexArray[k]));
           }
       }
   
   
   
   
       Compute_RotMat(InputData->RotMat_alpha, InputData->RotMat_beta, InputData->RotMat);
       LogSystem::write("Rotation Matrix: ");
       LogSystem::write(std::to_string(InputData->RotMat[0]) + " " + std::to_string(InputData->RotMat[3]) + " " + std::to_string(InputData->RotMat[6]));
       LogSystem::write(std::to_string(InputData->RotMat[1]) + " " + std::to_string(InputData->RotMat[4]) + " " + std::to_string(InputData->RotMat[7]));
       LogSystem::write(std::to_string(InputData->RotMat[2]) + " " + std::to_string(InputData->RotMat[5]) + " " + std::to_string(InputData->RotMat[8]));
       LogSystem::write("");
       LogSystem::write("");
       LogSystem::write("");
   
   
   
   
       // Check Polarization
       float P_norm = sqrt(pow(InputData->Polarization[0], 2) \
                         + pow(InputData->Polarization[1], 2) \
                         + pow(InputData->Polarization[2], 2));
       if(P_norm == 0){
           LogSystem::write("Error: Polarization magnitude is equal to zero!!");
       }
       if(P_norm != 0 && P_norm != 1){
           InputData->Polarization[0] = (InputData->Polarization[0])/P_norm;
           InputData->Polarization[1] = (InputData->Polarization[1])/P_norm;
           InputData->Polarization[2] = (InputData->Polarization[2])/P_norm;
           LogSystem::write("Polarization is automatically normalized to 1!");
       }
   
   
   
       LogSystem::write("##########################################################################################");
       LogSystem::write("## Stop - Input File Interpreter #########################################################");
       LogSystem::write("##########################################################################################");
       LogSystem::write("");
   
   
       bool ok = true;
   
       auto check_required = [&](auto& options)
       {
           for (auto& opt : options)
           {
               if (opt.required && !opt.found)
               {
                   LogSystem::write("Missing required option: " + opt.key);
                   ok = false;
               }
           }
       };
   
       check_required(Output_options);
       check_required(bool_options);
       check_required(int_options);
       check_required(float_options);
       check_required(string_options);
   
       if (!ok)
       {
           LogSystem::write(" ->-> Error in input file!");
       }
   
       return ok;
   
   }
   
   
