// File         : NuMagSANSlib_RotationDataObserver.h
// Author       : Michael Philipp ADAMS, M.Sc.
// Company      : University of Luxembourg
// Department   : Department of Physics and Materials Sciences
// Group        : NanoMagnetism Group
// Group Leader : Prof. Andreas Michels
// Version      : 25 May 2026
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

struct RotDataProperties{

    string GlobalFilePath;
    string GlobalFolderPath;

    int Number_Of_Elements;
    int Number_Of_Files;
    
};

// ###############################################################################################################################################
// helper functions ##############################################################################################################################
// ###############################################################################################################################################

void get_GlobalRotDataPath(std::string Local_RotDataPath,
                           RotDataProperties* RotDataProp){

    char tmp[PATH_MAX];
    getcwd(tmp, PATH_MAX);

    std::string tmp_string = tmp;

    // RotDataProp->GlobalFilePath = tmp_string + "/" + Local_RotDataPath;
    RotDataProp->GlobalFilePath = Local_RotDataPath;

    LogSystem::write("found Global RotDataPath: " + RotDataProp->GlobalFilePath);
}

void get_GlobalRotDataFolderPath(std::string Local_RotDataFolderPath,
                                 RotDataProperties* RotDataProp){

    char tmp[PATH_MAX];
    getcwd(tmp, PATH_MAX);

    std::string tmp_string = tmp;

    // RotDataProp->GlobalFolderPath = tmp_string + "/" + Local_RotDataFolderPath;
    RotDataProp->GlobalFolderPath = Local_RotDataFolderPath;

    LogSystem::write("found Global RotDataFolderPath: " + RotDataProp->GlobalFolderPath);
}

std::string RotDataLoopFilePath(RotDataProperties* RotDataProp,
                                int RotData_File_Index){

    return RotDataProp->GlobalFolderPath + "/RotData_" + std::to_string(RotData_File_Index) + ".csv";
}

bool FileExists(std::string filename){

    ifstream fin(filename);
    bool exists = fin.good();
    fin.close();
    return exists;
}

void SetActiveRotDataFile(RotDataProperties* RotDataProp,
                          int RotData_File_Index){

    RotDataProp->GlobalFilePath = RotDataLoopFilePath(RotDataProp, RotData_File_Index);
    LogSystem::write("active RotDataPath: " + RotDataProp->GlobalFilePath);
}


void NumberOfEntriesInRotationFile(int *NumberOfColumns,
                                   string filename){

    ifstream fin;
    fin.open(filename);

    std::string line;

    float alpha, beta, gamma;

    int line_counter = 0;
    int error_counter = 0;

    while(std::getline(fin, line)){
        std::istringstream ss(line);

        if(ss >> alpha >> beta >> gamma){
            line_counter += 1;
        } else {
            error_counter += 1;
        }
    }

    fin.close();

    *NumberOfColumns = line_counter;
}


// Routine that checks number of entries in RotationData file
bool RotData_Observer(std::string Local_RotDataPath,
                      RotDataProperties* RotDataProp){

    LogSystem::write("##########################################################################################");
    LogSystem::write("## Run - Rotation File Explorer ##########################################################");
    LogSystem::write("##########################################################################################");
    LogSystem::write("");

    bool CheckFlag = false;

    RotDataProp->Number_Of_Files = 1;

    // get global path of the RotationData file
    get_GlobalRotDataPath(Local_RotDataPath, RotDataProp);

    // count number of entries in the rotation data file
    NumberOfEntriesInRotationFile(&RotDataProp->Number_Of_Elements,
                                  RotDataProp->GlobalFilePath);

    LogSystem::write("Number of Entries: " + std::to_string(RotDataProp->Number_Of_Elements));
    LogSystem::write("");

    LogSystem::write("##########################################################################################");
    LogSystem::write("## Stop - Rotation File Explorer #########################################################");
    LogSystem::write("##########################################################################################");
    LogSystem::write("");
    LogSystem::write("");

    if(RotDataProp->Number_Of_Elements != 0){
        CheckFlag = true;
    }

    return CheckFlag;
}

bool RotDataLoop_Observer(std::string Local_RotDataFolderPath,
                          std::vector<int> RotDataLoop_IndexArray,
                          RotDataProperties* RotDataProp){

    LogSystem::write("##########################################################################################");
    LogSystem::write("## Run - Rotation Folder Explorer ########################################################");
    LogSystem::write("##########################################################################################");
    LogSystem::write("");

    bool CheckFlag = true;

    if(RotDataLoop_IndexArray.empty()){
        LogSystem::write("Error: RotDataLoop has no active RotData indices.");
        CheckFlag = false;
    }

    if(Local_RotDataFolderPath == ""){
        LogSystem::write("Error: RotDataPath is empty while RotDataLoop is active.");
        CheckFlag = false;
    }

    get_GlobalRotDataFolderPath(Local_RotDataFolderPath, RotDataProp);

    RotDataProp->Number_Of_Files = RotDataLoop_IndexArray.size();
    RotDataProp->Number_Of_Elements = 0;

    for(int k = 0; k < RotDataLoop_IndexArray.size(); k++){
        int RotData_File_Index = RotDataLoop_IndexArray[k];
        std::string filename = RotDataLoopFilePath(RotDataProp, RotData_File_Index);

        if(!FileExists(filename)){
            LogSystem::write("Error: missing RotData file: " + filename);
            CheckFlag = false;
            continue;
        }

        int Number_Of_Elements_tmp = 0;
        NumberOfEntriesInRotationFile(&Number_Of_Elements_tmp, filename);

        LogSystem::write("RotData_" + std::to_string(RotData_File_Index) + ".csv entries: "
                         + std::to_string(Number_Of_Elements_tmp));

        if(Number_Of_Elements_tmp == 0){
            LogSystem::write("Error: RotData file contains zero entries: " + filename);
            CheckFlag = false;
        }

        if(RotDataProp->Number_Of_Elements == 0){
            RotDataProp->Number_Of_Elements = Number_Of_Elements_tmp;
        }else if(RotDataProp->Number_Of_Elements != Number_Of_Elements_tmp){
            LogSystem::write("Error: RotData files do not contain the same number of entries.");
            CheckFlag = false;
        }
    }

    if(!RotDataLoop_IndexArray.empty()){
        SetActiveRotDataFile(RotDataProp, RotDataLoop_IndexArray[0]);
    }

    LogSystem::write("Number of RotData files: " + std::to_string(RotDataProp->Number_Of_Files));
    LogSystem::write("Number of Entries: " + std::to_string(RotDataProp->Number_Of_Elements));
    LogSystem::write("");

    LogSystem::write("##########################################################################################");
    LogSystem::write("## Stop - Rotation Folder Explorer #######################################################");
    LogSystem::write("##########################################################################################");
    LogSystem::write("");
    LogSystem::write("");

    return CheckFlag;
}
