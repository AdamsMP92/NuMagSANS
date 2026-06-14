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
    bool NuclearMagnetic = false;
    bool SpinFlip = false;
    bool Chiral = false;
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

/**
 * \ingroup InputParameters
 * \struct InputFileData
 * \brief Central configuration container for a NuMagSANS simulation run.
 *
 * This structure stores all parameters parsed from the user input file,
 * including geometry, discretization, physical constants, polarization
 * settings, and output control flags.
 */
struct InputFileData {

    /// Path to nuclear real-space data directory
    string NucDataPath;

    /// Path to magnetic real-space data directory
    string MagDataPath;

    /// CSV file containing structural information
    string StructDataFilename;

    /// CSV file containing rotational information
    string RotDataFilename;

    /// Directory containing StructData_1.csv, StructData_2.csv, ... files
    string StructDataPath = "";

    /// Directory containing RotData_1.csv, RotData_2.csv, ... files
    string RotDataPath = "";

    /// Skip repeated MagData dimension checks for faster data import
    bool FastLoad_flag = false;

    /// Replicate a single local nuclear object during import
    bool NucData_ReplicationImport_flag = false;

    /// Effective number of nuclear objects generated from the single template object
    int NucData_NumberOfReplications = 1;

    /// Replicate a single local magnetic object during import
    bool MagData_ReplicationImport_flag = false;

    /// Effective number of magnetic objects generated from the single template object
    int MagData_NumberOfReplications = 1;

    /// Output directory for computed SANS data
    string SANSDataFoldername;

    /// Fourier computation approach (e.g., "Direct", "FFT")
    string Fourier_Approach;

    /// Activate loop mode over dataset indices
    bool Loop_Modus;

    /// Activate inner loop mode over rotation data files
    bool RotDataLoop_flag = false;

    /// Activate inner loop mode over structure data files
    bool StructDataLoop_flag = false;

    /// Starting index for loop mode
    int Loop_From;

    /// Final index for loop mode
    int Loop_To;

    /// Starting index for rotation data loop mode
    int RotDataLoop_From = 1;

    /// Final index for rotation data loop mode
    int RotDataLoop_To = 1;

    /// Starting index for structure data loop mode
    int StructDataLoop_From = 1;

    /// Final index for structure data loop mode
    int StructDataLoop_To = 1;

    /// User-defined structure-data selection string
    string StructData_User_Selection;

    /// Parsed structure-data index list
    std::vector<int> StructDataLoop_IndexArray;

    /// User-defined rotation-data selection string
    string RotData_User_Selection;

    /// Parsed rotation-data index list
    std::vector<int> RotDataLoop_IndexArray;

    /// User-defined dataset selection string
    string User_Selection;

    /// Parsed dataset index list
    std::vector<int> User_Selection_IndexArray;

    // --- Discretization parameters ---

    /// Number of q-points
    int N_q;

    /// Number of azimuthal angles
    int N_theta;

    /// Number of radial points for correlation functions
    int N_r;

    /// Number of angular alpha points
    int N_alpha;

    /// Minimum scattering vector magnitude [1/nm]
    float q_min = 0.0;

    /// Maximum scattering vector magnitude [1/nm]
    float q_max;

    /// Maximum real-space radius [nm]
    float r_max;

    // --- Physical parameters ---

    /// Total scattering volume [nm^3]
    float Scattering_Volume_V;

    /// Nuclear scattering length density of a cell
    float cell_nuclear_sld;

    /// Saturation magnetization per cell
    float cell_magnetization;

    /// Cuboid cell size in x-direction [nm]
    float cuboid_cell_size_x;

    /// Cuboid cell size in y-direction [nm]
    float cuboid_cell_size_y;

    /// Cuboid cell size in z-direction [nm]
    float cuboid_cell_size_z;

    // --- Rotation matrix ---

    /// Euler rotation angle alpha [rad]
    float RotMat_alpha;

    /// Euler rotation angle beta [rad]
    float RotMat_beta;

    /// 3x3 rotation matrix (row-major)
    float RotMat[9];

    /// Unit conversion factor for coordinate scaling
    float XYZ_Unit_Factor;

    // --- Polarization vector ---

    /// Polarization vector (Px, Py, Pz)
    float Polarization[3];

    // --- Activation flags ---

    /// Enable nuclear data input
    bool NucData_activate_flag;

    /// Enable magnetic data input
    bool MagData_activate_flag;

    /// Enable structural data input
    bool StructData_activate_flag;

    /// Enable rotation data input
    bool RotData_activate_flag;

    /// Exclude zero-moment cells from computation
    bool ExcludeZeroMoments_flag;

    /// Validate input file before execution
    bool Check_InputFile_Flag;

    /// Maximum Fourier summation index
    int k_max;

    /// Activate angular spectrum evaluation
    bool AngularSpec_activate_flag;

    /// Output full Fourier correlation matrix
    bool output_fourier_correlation_matrix_flag;

    /*
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
    */
    OutFlagQuant OutFlags{};
};

inline bool any_active(const OutFlag& f) {
    return f.Nuclear || f.Unpolarized || f.Polarized || f.NuclearMagnetic || f.SpinFlip || f.Chiral || f.PM_SpinFlip ||
           f.MP_SpinFlip || f.PP_NonSpinFlip || f.MM_NonSpinFlip || f.P_SANSPOL || f.M_SANSPOL;
}

// zy-rotation matrix where alpha is the polar angle and beta is the azimuth angle
void Compute_RotMat(float alpha, float beta, float* RotMat) {

    // ordering of the rotation matrix:
    // RotMat = [RotMat[0], RotMat[3], RotMat[6]; ...
    //           RotMat[1], RotMat[4], RotMat[7]; ...
    //           RotMat[2], RotMat[5], RotMat[8]]

    alpha = alpha * M_PI / 180.0;
    beta = beta * M_PI / 180.0;

    RotMat[0] = cosf(alpha) * cosf(beta);
    RotMat[1] = cosf(alpha) * sinf(beta);
    RotMat[2] = -sinf(alpha);
    RotMat[3] = -sinf(beta);
    RotMat[4] = cosf(beta);
    RotMat[5] = 0;
    RotMat[6] = cosf(beta) * sinf(alpha);
    RotMat[7] = sinf(alpha) * sinf(beta);
    RotMat[8] = cosf(alpha);
}

static std::string trim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    size_t end = s.find_last_not_of(" \t\r\n");

    if (start == std::string::npos)
        return "";

    return s.substr(start, end - start + 1);
}

static bool try_extract_value(const std::string& line, const std::string& property, std::string& value) {
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
    while (pos < trimmed.size() && std::isspace(static_cast<unsigned char>(trimmed[pos])))
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
    value = trim(trimmed.substr(pos + 1, semi_pos - pos - 1));

    return true;
}

void parseStringNoCout(const std::string& line, const std::string& property, std::string& variable, bool& flag) {
    if (!try_extract_value(line, property, variable))
        return;

    flag = true;
}

static std::vector<int> parse_int_list(const std::string& s) {
    std::vector<int> result;

    std::string trimmed = trim(s);

    if (trimmed.front() != '{' || trimmed.back() != '}')
        throw std::runtime_error("Invalid integer list format");

    std::string content = trimmed.substr(1, trimmed.size() - 2);

    std::stringstream ss(content);
    std::string item;

    while (std::getline(ss, item, ',')) {
        result.push_back(std::stoi(trim(item)));
    }

    return result;
}

// Interpreting User_Selection to integer array ##########################################
void User_Selection_To_Int_Array(int* Idx, int Number_of_Comma, string my_string) {

    int Symbol_Idx[Number_of_Comma + 2];
    string Symbol_Start = "{";
    string Symbol_End = "}";
    string Symbol_Comma = ",";
    int Symbol_Counter = 0;
    for (int i = 0; i < my_string.size(); i++) {
        if (my_string.at(i) == Symbol_Start.at(0) || my_string.at(i) == Symbol_End.at(0) ||
            my_string.at(i) == Symbol_Comma.at(0)) {
            Symbol_Idx[Symbol_Counter] = i;
            Symbol_Counter += 1;
        }
    }

    string Sub_String;
    int space_symbol_index;
    for (int i = 0; i < Number_of_Comma + 1; i++) {
        Sub_String = my_string.substr(Symbol_Idx[i] + 1, Symbol_Idx[i + 1] - Symbol_Idx[i] - 1);
        while (Sub_String.find(" ") <= Sub_String.size()) {
            space_symbol_index = Sub_String.find(" ");
            Sub_String.erase(space_symbol_index, 1);
        }
        Idx[i] = stoi(Sub_String);
    }
}

// Find Number of Comma in string #########################################################
void Number_Of_Comma_In_String(int* N, string my_string) {
    int String_Length = my_string.size();
    int Comma_Counter = 0;
    string Comma = ",";

    for (int i = 0; i < String_Length; i++) {
        if (my_string.at(i) == Comma.at(0)) {
            Comma_Counter += 1;
        }
    }
    *N = Comma_Counter;
}

// Pre-Parsing for the initialization of the LogFile System
void InitializeLogSystem(string filename) {

    ifstream fin;
    fin.open(filename);
    string line;
    string string_to_be_interpreted;
    bool checkFlag = false;

    string SANSDataFoldername;

    while (getline(fin, line)) {
        parseStringNoCout(line, "foldernameSANSData", SANSDataFoldername, checkFlag);
    }
    fin.close();

    if (checkFlag) {
        mkdir((SANSDataFoldername + "/").c_str(), 0777);
        LogSystem::initLog(SANSDataFoldername + "/NuMagSANSlog.txt");
        // LogSystem::write("LogFile Initialized");
    } else {
        std::cerr << "Error: 'foldernameSANSData' not found in input file '" << filename
                  << "'. Logging not initialized.\n";
        std::cerr << "Aborting program.\n";
        std::exit(EXIT_FAILURE);
    }
}

// template for options that are parsed by the InputFileInterpreter
template <typename T> struct Option {
    std::string key;
    T* target;
    bool required;
    bool found = false;
};

bool ReadCSV__Input_File_Interpreter(string filename, InputFileData* InputData) {

    LogSystem::write("##########################################################################################");
    LogSystem::write("## Run - Input File Interpreter ##########################################################");
    LogSystem::write("##########################################################################################");

    ifstream fin(filename);
    if (!fin.is_open()) {
        LogSystem::write("Could not open input file.");
        return false;
    }

    // generalisation of the output options
    std::vector<PrefixEntry> prefixes = {{"Nuclear", &OutFlag::Nuclear},
                                         {"Unpolarized", &OutFlag::Unpolarized},
                                         {"Polarized", &OutFlag::Polarized},
                                         {"NuclearMagnetic", &OutFlag::NuclearMagnetic},
                                         {"SpinFlip", &OutFlag::SpinFlip},
                                         {"Chiral", &OutFlag::Chiral},
                                         {"PM_SpinFlip", &OutFlag::PM_SpinFlip},
                                         {"MP_SpinFlip", &OutFlag::MP_SpinFlip},
                                         {"PP_NonSpinFlip", &OutFlag::PP_NonSpinFlip},
                                         {"MM_NonSpinFlip", &OutFlag::MM_NonSpinFlip},
                                         {"P_SANSPOL", &OutFlag::P_SANSPOL},
                                         {"M_SANSPOL", &OutFlag::M_SANSPOL}};

    std::vector<SuffixEntry> suffixes = {{"_2D", &OutFlagQuant::SANS2D},
                                         {"_1D", &OutFlagQuant::SANS1D},
                                         {"_Corr_2D", &OutFlagQuant::Corr2D},
                                         {"_Corr_1D", &OutFlagQuant::Corr1D},
                                         {"_PairDist_1D", &OutFlagQuant::PairDist1D}};

    std::vector<Option<bool>> Output_options;
    for (auto& p : prefixes) {
        for (auto& s : suffixes) {
            std::string key = std::string(p.name) + s.name;
            Output_options.push_back({key, &((InputData->OutFlags.*(s.quantity)).*(p.member)), true});
        }
    }

    // bool options //////////////////////////////
    std::vector<Option<bool>> bool_options = {

        {"Loop_Modus", &InputData->Loop_Modus, true},
        {"NucData_activate", &InputData->NucData_activate_flag, true},
        {"MagData_activate", &InputData->MagData_activate_flag, true},
        {"StructData_activate", &InputData->StructData_activate_flag, true},
        {"RotData_activate", &InputData->RotData_activate_flag, true},
        {"StructDataLoop", &InputData->StructDataLoop_flag, false},
        {"RotDataLoop", &InputData->RotDataLoop_flag, false},
        {"FastLoad", &InputData->FastLoad_flag, false},
        {"NucData_ReplicationImport", &InputData->NucData_ReplicationImport_flag, false},
        {"MagData_ReplicationImport", &InputData->MagData_ReplicationImport_flag, false},
        {"Exclude_Zero_Moments", &InputData->ExcludeZeroMoments_flag, true},
        {"Angular_Spec", &InputData->AngularSpec_activate_flag, true},
        {"Fourier_Gamma", &InputData->output_fourier_correlation_matrix_flag, true}

    };

    // int options /////////////////////////////////
    std::vector<Option<int>> int_options = {

        {"Loop_From", &InputData->Loop_From, true},
        {"Loop_To", &InputData->Loop_To, true},
        {"StructDataLoop_From", &InputData->StructDataLoop_From, false},
        {"StructDataLoop_To", &InputData->StructDataLoop_To, false},
        {"RotDataLoop_From", &InputData->RotDataLoop_From, false},
        {"RotDataLoop_To", &InputData->RotDataLoop_To, false},
        {"NucData_NumberOfReplications", &InputData->NucData_NumberOfReplications, false},
        {"MagData_NumberOfReplications", &InputData->MagData_NumberOfReplications, false},
        {"Number_Of_q_Points", &InputData->N_q, true},
        {"Number_Of_theta_Points", &InputData->N_theta, true},
        {"Number_Of_r_Points", &InputData->N_r, true},
        {"Number_Of_alpha_Points", &InputData->N_alpha, true},
        {"k_max", &InputData->k_max, true}

    };

    // float options ///////////////////////////////
    std::vector<Option<float>> float_options = {

        {"q_min", &InputData->q_min, false},
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
        {"RotDataFilename", &InputData->RotDataFilename, true},
        {"StructDataPath", &InputData->StructDataPath, false},
        {"RotDataPath", &InputData->RotDataPath, false},
        {"foldernameSANSData", &InputData->SANSDataFoldername, true},
        {"Fourier_Approach", &InputData->Fourier_Approach, true},
        {"StructData_User_Selection", &InputData->StructData_User_Selection, false},
        {"RotData_User_Selection", &InputData->RotData_User_Selection, false},
        {"User_Selection", &InputData->User_Selection, true}

    };

    std::string line;
    while (getline(fin, line)) {
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
        if (opt.key == "User_Selection" && opt.found) {
            InputData->User_Selection_IndexArray = parse_int_list(InputData->User_Selection);

            LogSystem::write("Check UserSelection entries that are transferred to integer array:");

            for (size_t k = 0; k < InputData->User_Selection_IndexArray.size(); ++k) {
                LogSystem::write(" UserSelection: " + std::to_string(k) + " : " +
                                 std::to_string(InputData->User_Selection_IndexArray[k]));
            }
        }

    for (auto& opt : string_options)
        if (opt.key == "StructData_User_Selection" && opt.found) {
            InputData->StructDataLoop_IndexArray = parse_int_list(InputData->StructData_User_Selection);

            LogSystem::write("Check StructDataUserSelection entries that are transferred to integer array:");

            for (size_t k = 0; k < InputData->StructDataLoop_IndexArray.size(); ++k) {
                LogSystem::write(" StructDataUserSelection: " + std::to_string(k) + " : " +
                                 std::to_string(InputData->StructDataLoop_IndexArray[k]));
            }
        }

    if (InputData->StructDataLoop_flag && InputData->StructDataLoop_IndexArray.empty()) {
        for (int k = InputData->StructDataLoop_From; k <= InputData->StructDataLoop_To; k++) {
            InputData->StructDataLoop_IndexArray.push_back(k);
        }
    }

    for (auto& opt : string_options)
        if (opt.key == "RotData_User_Selection" && opt.found) {
            InputData->RotDataLoop_IndexArray = parse_int_list(InputData->RotData_User_Selection);

            LogSystem::write("Check RotDataUserSelection entries that are transferred to integer array:");

            for (size_t k = 0; k < InputData->RotDataLoop_IndexArray.size(); ++k) {
                LogSystem::write(" RotDataUserSelection: " + std::to_string(k) + " : " +
                                 std::to_string(InputData->RotDataLoop_IndexArray[k]));
            }
        }

    if (InputData->RotDataLoop_flag && InputData->RotDataLoop_IndexArray.empty()) {
        for (int k = InputData->RotDataLoop_From; k <= InputData->RotDataLoop_To; k++) {
            InputData->RotDataLoop_IndexArray.push_back(k);
        }
    }

    Compute_RotMat(InputData->RotMat_alpha, InputData->RotMat_beta, InputData->RotMat);
    LogSystem::write("Rotation Matrix: ");
    LogSystem::write(std::to_string(InputData->RotMat[0]) + " " + std::to_string(InputData->RotMat[3]) + " " +
                     std::to_string(InputData->RotMat[6]));
    LogSystem::write(std::to_string(InputData->RotMat[1]) + " " + std::to_string(InputData->RotMat[4]) + " " +
                     std::to_string(InputData->RotMat[7]));
    LogSystem::write(std::to_string(InputData->RotMat[2]) + " " + std::to_string(InputData->RotMat[5]) + " " +
                     std::to_string(InputData->RotMat[8]));
    LogSystem::write("");
    LogSystem::write("");
    LogSystem::write("");

    // Check Polarization
    float P_norm = sqrtf(powf(InputData->Polarization[0], 2) + powf(InputData->Polarization[1], 2) +
                         powf(InputData->Polarization[2], 2));
    if (P_norm == 0) {
        LogSystem::write("Error: Polarization magnitude is equal to zero!!");
    }
    if (P_norm != 0 && P_norm != 1) {
        InputData->Polarization[0] = (InputData->Polarization[0]) / P_norm;
        InputData->Polarization[1] = (InputData->Polarization[1]) / P_norm;
        InputData->Polarization[2] = (InputData->Polarization[2]) / P_norm;
        LogSystem::write("Polarization is automatically normalized to 1!");
    }

    bool ReplicationImport_CheckFlag = true;
    if (InputData->MagData_ReplicationImport_flag && !InputData->MagData_activate_flag) {
        LogSystem::write("Error: MagData_ReplicationImport requires MagData_activate = 1.");
        ReplicationImport_CheckFlag = false;
    }
    if (InputData->NucData_ReplicationImport_flag && !InputData->NucData_activate_flag) {
        LogSystem::write("Error: NucData_ReplicationImport requires NucData_activate = 1.");
        ReplicationImport_CheckFlag = false;
    }
    if (InputData->MagData_ReplicationImport_flag && InputData->MagData_NumberOfReplications <= 0) {
        LogSystem::write("Error: MagData_NumberOfReplications must be larger than zero.");
        ReplicationImport_CheckFlag = false;
    }
    if (InputData->NucData_ReplicationImport_flag && InputData->NucData_NumberOfReplications <= 0) {
        LogSystem::write("Error: NucData_NumberOfReplications must be larger than zero.");
        ReplicationImport_CheckFlag = false;
    }

    bool ScatteringGrid_CheckFlag = true;
    if (InputData->N_q <= 1) {
        LogSystem::write("Error: Number_Of_q_Points must be larger than one.");
        ScatteringGrid_CheckFlag = false;
    }
    if (InputData->q_min < 0.0) {
        LogSystem::write("Error: q_min must be larger than or equal to zero.");
        ScatteringGrid_CheckFlag = false;
    }
    if (InputData->q_max <= InputData->q_min) {
        LogSystem::write("Error: q_max must be larger than q_min.");
        ScatteringGrid_CheckFlag = false;
    }

    LogSystem::write("##########################################################################################");
    LogSystem::write("## Stop - Input File Interpreter #########################################################");
    LogSystem::write("##########################################################################################");
    LogSystem::write("");

    bool ok = true;

    auto check_required = [&](auto& options) {
        for (auto& opt : options) {
            if (opt.required && !opt.found) {
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

    if (!ok) {
        LogSystem::write(" ->-> Error in input file!");
    }

    return ok && ReplicationImport_CheckFlag && ScatteringGrid_CheckFlag;
}
