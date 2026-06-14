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

struct StructDataProperties {

    string GlobalFilePath;
    string GlobalFolderPath;

    int Number_Of_Elements;
    int Number_Of_Files;
};

// ###############################################################################################################################################
// helper functions
// ##############################################################################################################################
// ###############################################################################################################################################

void get_GlobalStructDataPath(std::string Local_StructDataPath, StructDataProperties* StructDataProp) {

    char tmp[PATH_MAX];
    getcwd(tmp, PATH_MAX); // Get the current working directory
    std::string tmp_string = tmp;
    // StructDataProp->GlobalFilePath = tmp_string + "/" + Local_StructDataPath;
    StructDataProp->GlobalFilePath = Local_StructDataPath;
    LogSystem::write("found Global StructDataPath: " + StructDataProp->GlobalFilePath);
}

void get_GlobalStructDataFolderPath(std::string Local_StructDataFolderPath, StructDataProperties* StructDataProp) {

    char tmp[PATH_MAX];
    getcwd(tmp, PATH_MAX);

    std::string tmp_string = tmp;

    // StructDataProp->GlobalFolderPath = tmp_string + "/" + Local_StructDataFolderPath;
    StructDataProp->GlobalFolderPath = Local_StructDataFolderPath;
    LogSystem::write("found Global StructDataFolderPath: " + StructDataProp->GlobalFolderPath);
}

std::string StructDataLoopFilePath(StructDataProperties* StructDataProp, int StructData_File_Index) {

    return StructDataProp->GlobalFolderPath + "/StructData_" + std::to_string(StructData_File_Index) + ".csv";
}

bool StructFileExists(std::string filename) {

    ifstream fin(filename);
    bool exists = fin.good();
    fin.close();
    return exists;
}

void SetActiveStructDataFile(StructDataProperties* StructDataProp, int StructData_File_Index) {

    StructDataProp->GlobalFilePath = StructDataLoopFilePath(StructDataProp, StructData_File_Index);
    LogSystem::write("active StructDataPath: " + StructDataProp->GlobalFilePath);
}

void NumberOfEntriesInStructureFile(int* NumberOfColumns, string filename) {

    ifstream fin;
    fin.open(filename);
    std::string line;

    float x, y, z;

    int line_counter = 0;
    int error_counter = 0;

    while (std::getline(fin, line)) {
        std::istringstream ss(line);
        if (ss >> x >> y >> z) {
            line_counter += 1;
        } else {
            error_counter++;
        }
    }
    fin.close();
    *NumberOfColumns = line_counter;
}

// Routine that checks number of subfolders in MagData directory
bool StructData_Observer(std::string Local_StructDataPath, StructDataProperties* StructDataProp) {

    LogSystem::write("##########################################################################################");
    LogSystem::write("## Run - Structure File Explorer #########################################################");
    LogSystem::write("##########################################################################################");
    LogSystem::write("");

    bool CheckFlag = false;

    StructDataProp->Number_Of_Files = 1;

    // get global path of the MagData directory
    get_GlobalStructDataPath(Local_StructDataPath, StructDataProp);

    // count number of entries in the structure data file
    NumberOfEntriesInStructureFile(&StructDataProp->Number_Of_Elements, StructDataProp->GlobalFilePath);
    LogSystem::write("Number of Entries: " + std::to_string(StructDataProp->Number_Of_Elements));
    LogSystem::write("");

    LogSystem::write("##########################################################################################");
    LogSystem::write("## Stop - Structure File Explorer ########################################################");
    LogSystem::write("##########################################################################################");
    LogSystem::write("");
    LogSystem::write("");

    if (StructDataProp->Number_Of_Elements != 0) {
        CheckFlag = true;
    }

    return CheckFlag;
}

bool StructDataLoop_Observer(std::string Local_StructDataFolderPath, std::vector<int> StructDataLoop_IndexArray,
                             StructDataProperties* StructDataProp) {

    LogSystem::write("##########################################################################################");
    LogSystem::write("## Run - Structure Folder Explorer #######################################################");
    LogSystem::write("##########################################################################################");
    LogSystem::write("");

    bool CheckFlag = true;

    if (StructDataLoop_IndexArray.empty()) {
        LogSystem::write("Error: StructDataLoop has no active StructData indices.");
        CheckFlag = false;
    }

    if (Local_StructDataFolderPath == "") {
        LogSystem::write("Error: StructDataPath is empty while StructDataLoop is active.");
        CheckFlag = false;
    }

    get_GlobalStructDataFolderPath(Local_StructDataFolderPath, StructDataProp);

    StructDataProp->Number_Of_Files = StructDataLoop_IndexArray.size();
    StructDataProp->Number_Of_Elements = 0;

    for (int k = 0; k < StructDataLoop_IndexArray.size(); k++) {
        int StructData_File_Index = StructDataLoop_IndexArray[k];
        std::string filename = StructDataLoopFilePath(StructDataProp, StructData_File_Index);

        if (!StructFileExists(filename)) {
            LogSystem::write("Error: missing StructData file: " + filename);
            CheckFlag = false;
            continue;
        }

        int Number_Of_Elements_tmp = 0;
        NumberOfEntriesInStructureFile(&Number_Of_Elements_tmp, filename);

        LogSystem::write("StructData_" + std::to_string(StructData_File_Index) +
                         ".csv entries: " + std::to_string(Number_Of_Elements_tmp));

        if (Number_Of_Elements_tmp == 0) {
            LogSystem::write("Error: StructData file contains zero entries: " + filename);
            CheckFlag = false;
        }

        if (StructDataProp->Number_Of_Elements == 0) {
            StructDataProp->Number_Of_Elements = Number_Of_Elements_tmp;
        } else if (StructDataProp->Number_Of_Elements != Number_Of_Elements_tmp) {
            LogSystem::write("Error: StructData files do not contain the same number of entries.");
            CheckFlag = false;
        }
    }

    if (!StructDataLoop_IndexArray.empty()) {
        SetActiveStructDataFile(StructDataProp, StructDataLoop_IndexArray[0]);
    }

    LogSystem::write("Number of StructData files: " + std::to_string(StructDataProp->Number_Of_Files));
    LogSystem::write("Number of Entries: " + std::to_string(StructDataProp->Number_Of_Elements));
    LogSystem::write("");

    LogSystem::write("##########################################################################################");
    LogSystem::write("## Stop - Structure Folder Explorer ######################################################");
    LogSystem::write("##########################################################################################");
    LogSystem::write("");
    LogSystem::write("");

    return CheckFlag;
}
