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
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>

using namespace std;

// ===============================================================
// DATA STRUCT
// ===============================================================

struct InputFileData {

    // folders
    string NucDataPath;
    string MagDataPath;
    string StructDataFilename;
    string SANSDataFoldername;
    string Fourier_Approach;

    // loop
    bool Loop_Modus = false;
    int Loop_From = 0;
    int Loop_To = 0;

    // selection
    string User_Selection;
    vector<int> User_Selection_IndexArray;

    // grid
    int N_q = 0;
    int N_theta = 0;
    int N_r = 0;
    int N_alpha = 0;
    float q_max = 0;
    float r_max = 0;

    // micromagnetics
    float Scattering_Volume_V = 0;
    float cell_nuclear_sld = 0;
    float cell_magnetization = 0;
    float cuboid_cell_size_x = 0;
    float cuboid_cell_size_y = 0;
    float cuboid_cell_size_z = 0;

    // rotation
    float RotMat_alpha = 0;
    float RotMat_beta = 0;
    float RotMat[9] = {0};

    float XYZ_Unit_Factor = 1.0f;
    float Polarization[3] = {0,0,1};

    // flags
    bool NucData_activate_flag = false;
    bool MagData_activate_flag = false;
    bool StructData_activate_flag = false;
    bool ExcludeZeroMoments_flag = false;

    int k_max = 0;
    bool AngularSpec_activate_flag = false;

    // outputs 
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

// ===============================================================
// UTILITIES
// ===============================================================

inline void trim(string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    size_t end   = s.find_last_not_of(" \t\r\n");

    if (start == string::npos)
        s = "";
    else
        s = s.substr(start, end - start + 1);
}

inline bool extractKeyValue(const string& line,
                            string& key,
                            string& value)
{
    size_t eq = line.find('=');
    size_t sc = line.find(';');

    if (eq == string::npos || sc == string::npos || sc < eq)
        return false;

    key = line.substr(0, eq);
    value = line.substr(eq + 1, sc - eq - 1);

    trim(key);
    trim(value);

    return true;
}

// ===============================================================
// ROTATION MATRIX
// ===============================================================

inline void Compute_RotMat(float alpha, float beta, float* R) {

    alpha *= M_PI/180.0f;
    beta  *= M_PI/180.0f;

    R[0] = cos(alpha) * cos(beta);
    R[1] = cos(alpha) * sin(beta);
    R[2] = -sin(alpha);
    R[3] = -sin(beta);
    R[4] = cos(beta);
    R[5] = 0;
    R[6] = cos(beta) * sin(alpha);
    R[7] = sin(alpha) * sin(beta);
    R[8] = cos(alpha);
}

// ===============================================================
// MAIN PARSER
// ===============================================================

inline bool ReadCSV__Input_File_Interpreter(
    const string& filename,
    InputFileData& data)
{
    ifstream fin(filename);
    if (!fin)
        throw runtime_error("Cannot open input file.");

    LogSystem::write("=================================================");
    LogSystem::write("InputFileInterpreter2 started");
    LogSystem::write("=================================================");

    unordered_map<string, function<void(const string&)>> registry;
    unordered_set<string> seen_keys;

    // ----------------------------
    // Binder Helpers
    // ----------------------------

    auto bindFloat = [&](const string& key, float& var){
        registry[key] = [&](const string& v){
            var = stof(v);
            LogSystem::write(" -> " + key + " : " + to_string(var));
        };
    };

    auto bindInt = [&](const string& key, int& var){
        registry[key] = [&](const string& v){
            var = stoi(v);
            LogSystem::write(" -> " + key + " : " + to_string(var));
        };
    };

    auto bindBool = [&](const string& key, bool& var){
        registry[key] = [&](const string& v){
            var = static_cast<bool>(stoi(v));
            LogSystem::write(" -> " + key + " : " + string(var ? "true" : "false"));
        };
    };

    auto bindString = [&](const string& key, string& var){
        registry[key] = [&](const string& v){
            var = v;
            LogSystem::write(" -> " + key + " : " + var);
        };
    };

    // ----------------------------
    // REGISTRY
    // ----------------------------

    // folders
    bindString("NucDataPath", data.NucDataPath);
    bindString("MagDataPath", data.MagDataPath);
    bindString("StructDataFilename", data.StructDataFilename);
    bindString("foldernameSANSData", data.SANSDataFoldername);
    bindString("Fourier_Approach", data.Fourier_Approach);

    // loop
    bindBool("Loop_Modus", data.Loop_Modus);
    bindInt("Loop_From", data.Loop_From);
    bindInt("Loop_To", data.Loop_To);

    // selection
    registry["User_Selection"] = [&](const string& v){
        data.User_Selection = v;
        LogSystem::write(" -> User_Selection : " + v);
    };

    // grid
    bindInt("Number_Of_q_Points", data.N_q);
    bindInt("Number_Of_theta_Points", data.N_theta);
    bindInt("Number_Of_r_Points", data.N_r);
    bindInt("Number_Of_alpha_Points", data.N_alpha);
    bindFloat("q_max", data.q_max);
    bindFloat("r_max", data.r_max);

    // micromagnetics
    bindFloat("Scattering_Volume_V", data.Scattering_Volume_V);
    bindFloat("Cell_Nuclear_SLD", data.cell_nuclear_sld);
    bindFloat("Cell_Magnetization", data.cell_magnetization);
    bindFloat("Cuboid_Cell_Size_x", data.cuboid_cell_size_x);
    bindFloat("Cuboid_Cell_Size_y", data.cuboid_cell_size_y);
    bindFloat("Cuboid_Cell_Size_z", data.cuboid_cell_size_z);

    // rotation
    bindFloat("RotMat_alpha", data.RotMat_alpha);
    bindFloat("RotMat_beta", data.RotMat_beta);

    // misc
    bindFloat("XYZ_Unit_Factor", data.XYZ_Unit_Factor);
    bindFloat("Polarization_x", data.Polarization[0]);
    bindFloat("Polarization_y", data.Polarization[1]);
    bindFloat("Polarization_z", data.Polarization[2]);

    bindBool("NucData_activate", data.NucData_activate_flag);
    bindBool("MagData_activate", data.MagData_activate_flag);
    bindBool("StructData_activate", data.StructData_activate_flag);
    bindBool("Exclude_Zero_Moments", data.ExcludeZeroMoments_flag);

    bindInt("k_max", data.k_max);
    bindBool("Angular_Spec", data.AngularSpec_activate_flag);

    // ----------------------------
    // ALL OUTPUT FLAGS (vollständig)
    // ----------------------------

    bindBool("Fourier_Gamma", data.output_fourier_correlation_matrix_flag);

    // 2D
    bindBool("Nuclear_2D", data.output_unpolarized_nuclear_SANS_cross_section_2D_flag);
    bindBool("Unpolarized_2D", data.output_unpolarized_magnetic_SANS_cross_section_2D_flag);
    bindBool("Polarized_2D", data.output_polarized_magnetic_SANS_cross_section_2D_flag);
    bindBool("NuclearMagnetic_2D", data.output_nuclear_magnetic_SANS_cross_section_2D_flag);
    bindBool("AVG_SpinFlip_2D", data.output_spin_flip_magnetic_SANS_cross_section_2D_flag);
    bindBool("Chiral_2D", data.output_chiral_magnetic_SANS_cross_section_2D_flag);
    bindBool("PM_SpinFlip_2D", data.output_pm_spin_flip_magnetic_SANS_cross_section_2D_flag);
    bindBool("MP_SpinFlip_2D", data.output_mp_spin_flip_magnetic_SANS_cross_section_2D_flag);
    bindBool("PP_NonSpinFlip_2D", data.output_pp_non_spin_flip_magnetic_SANS_cross_section_2D_flag);
    bindBool("MM_NonSpinFlip_2D", data.output_mm_non_spin_flip_magnetic_SANS_cross_section_2D_flag);
    bindBool("P_SANSPOL_2D", data.output_p_sanspol_magnetic_SANS_cross_section_2D_flag);
    bindBool("M_SANSPOL_2D", data.output_m_sanspol_magnetic_SANS_cross_section_2D_flag);

    // 1D
    bindBool("Nuclear_1D", data.output_unpolarized_nuclear_SANS_cross_section_1D_flag);
    bindBool("Unpolarized_1D", data.output_unpolarized_magnetic_SANS_cross_section_1D_flag);
    bindBool("Polarized_1D", data.output_polarized_magnetic_SANS_cross_section_1D_flag);
    bindBool("NuclearMagnetic_1D", data.output_nuclear_magnetic_SANS_cross_section_1D_flag);
    bindBool("AVG_SpinFlip_1D", data.output_spin_flip_magnetic_SANS_cross_section_1D_flag);
    bindBool("Chiral_1D", data.output_chiral_magnetic_SANS_cross_section_1D_flag);
    bindBool("PM_SpinFlip_1D", data.output_pm_spin_flip_SANS_cross_section_1D_flag);
    bindBool("MP_SpinFlip_1D", data.output_mp_spin_flip_SANS_cross_section_1D_flag);
    bindBool("PP_NonSpinFlip_1D", data.output_pp_non_spin_flip_SANS_cross_section_1D_flag);
    bindBool("MM_NonSpinFlip_1D", data.output_mm_non_spin_flip_SANS_cross_section_1D_flag);
    bindBool("P_SANSPOL_1D", data.output_p_sanspol_cross_section_1D_flag);
    bindBool("M_SANSPOL_1D", data.output_m_sanspol_cross_section_1D_flag);

    // PairDist 1D
    bindBool("Nuclear_PairDist_1D", data.output_nuclear_pair_distance_distribution_1D_flag);
    bindBool("Unpolarized_PairDist_1D", data.output_unpolarized_pair_distance_distribution_1D_flag);
    bindBool("Polarized_PairDist_1D", data.output_polarized_pair_distance_distribution_1D_flag);
    bindBool("NucMag_PairDist_1D", data.output_nuclear_magnetic_pair_distance_distribution_1D_flag);
    bindBool("SpinFlip_PairDist_1D", data.output_spin_flip_pair_distance_distribution_1D_flag);
    bindBool("Chiral_PairDist_1D", data.output_chiral_pair_distance_distribution_1D_flag);
    bindBool("PM_SpinFlip_PairDist_1D", data.output_pm_spin_flip_pair_distance_distribution_1D_flag);
    bindBool("MP_SpinFlip_PairDist_1D", data.output_mp_spin_flip_pair_distance_distribution_1D_flag);
    bindBool("PP_NonSpinFlip_PairDist_1D", data.output_pp_non_spin_flip_pair_distance_distribution_1D_flag);
    bindBool("MM_NonSpinFlip_PairDist_1D", data.output_mm_non_spin_flip_pair_distance_distribution_1D_flag);
    bindBool("P_SANSPOL_PairDist_1D", data.output_p_sanspol_pair_distance_distribution_1D_flag);
    bindBool("M_SANSPOL_PairDist_1D", data.output_m_sanspol_pair_distance_distribution_1D_flag);

    // Corr 1D
    bindBool("Nuclear_Corr_1D", data.output_nuclear_correlation_function_1D_flag);
    bindBool("Unpolarized_Corr_1D", data.output_unpolarized_correlation_function_1D_flag);
    bindBool("Polarized_Corr_1D", data.output_polarized_correlation_function_1D_flag);
    bindBool("NucMag_Corr_1D", data.output_nuclear_magnetic_correlation_function_1D_flag);
    bindBool("SpinFlip_Corr_1D", data.output_spin_flip_correlation_function_1D_flag);
    bindBool("Chiral_Corr_1D", data.output_chiral_correlation_function_1D_flag);
    bindBool("PM_SpinFlip_Corr_1D", data.output_pm_spin_flip_correlation_function_1D_flag);
    bindBool("MP_SpinFlip_Corr_1D", data.output_mp_spin_flip_correlation_function_1D_flag);
    bindBool("PP_NonSpinFlip_Corr_1D", data.output_pp_non_spin_flip_correlation_function_1D_flag);
    bindBool("MM_NonSpinFlip_Corr_1D", data.output_mm_non_spin_flip_correlation_function_1D_flag);
    bindBool("P_SANSPOL_Corr_1D", data.output_p_sanspol_correlation_function_1D_flag);
    bindBool("M_SANSPOL_Corr_1D", data.output_m_sanspol_correlation_function_1D_flag);

    // Corr 2D
    bindBool("Nuclear_Corr_2D", data.output_nuclear_correlation_function_2D_flag);
    bindBool("Unpolarized_Corr_2D", data.output_unpolarized_correlation_function_2D_flag);
    bindBool("Polarized_Corr_2D", data.output_polarized_correlation_function_2D_flag);
    bindBool("NuclearMagnetic_Corr_2D", data.output_nuclear_magnetic_correlation_function_2D_flag);
    bindBool("SpinFlip_Corr_2D", data.output_spin_flip_correlation_function_2D_flag);
    bindBool("Chiral_Corr_2D", data.output_chiral_correlation_function_2D_flag);
    bindBool("PM_SpinFlip_Corr_2D", data.output_pm_spin_flip_correlation_function_2D_flag);
    bindBool("MP_SpinFlip_Corr_2D", data.output_mp_spin_flip_correlation_function_2D_flag);
    bindBool("PP_NonSpinFlip_Corr_2D", data.output_pp_non_spin_flip_correlation_function_2D_flag);
    bindBool("MM_NonSpinFlip_Corr_2D", data.output_mm_non_spin_flip_correlation_function_2D_flag);
    bindBool("P_SANSPOL_Corr_2D", data.output_p_sanspol_correlation_function_2D_flag);
    bindBool("M_SANSPOL_Corr_2D", data.output_m_sanspol_correlation_function_2D_flag);

    // ----------------------------
    // PARSING LOOP
    // ----------------------------

    string line, key, value;

    while (getline(fin, line)) {

        if (!extractKeyValue(line, key, value))
            continue;

        auto it = registry.find(key);

        if (it != registry.end()) {
            it->second(value);
            seen_keys.insert(key);
        }
        else {
            LogSystem::write("Warning: Unknown key -> " + key);
        }
    }

    fin.close();

    // ----------------------------
    // User_Selection parsing
    // ----------------------------

    if (!data.User_Selection.empty()) {
        data.User_Selection_IndexArray.clear();

        string s = data.User_Selection;
        s.erase(remove(s.begin(), s.end(), '{'), s.end());
        s.erase(remove(s.begin(), s.end(), '}'), s.end());

        stringstream ss(s);
        string item;

        while (getline(ss, item, ',')) {
            trim(item);
            if (!item.empty())
                data.User_Selection_IndexArray.push_back(stoi(item));
        }
    }

    // ----------------------------
    // Post Processing
    // ----------------------------

    Compute_RotMat(data.RotMat_alpha,
                   data.RotMat_beta,
                   data.RotMat);

    float norm = sqrt(
        data.Polarization[0]*data.Polarization[0] +
        data.Polarization[1]*data.Polarization[1] +
        data.Polarization[2]*data.Polarization[2]);

    if (norm == 0)
        throw runtime_error("Polarization vector is zero.");

    if (fabs(norm - 1.0f) > 1e-6f) {
        data.Polarization[0] /= norm;
        data.Polarization[1] /= norm;
        data.Polarization[2] /= norm;
        LogSystem::write("Polarization automatically normalized.");
    }

    LogSystem::write("InputFileInterpreter2 finished.");
    LogSystem::write("=================================================");

    return true;
}