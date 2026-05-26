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

    int Number_Of_Elements;
    
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

    cout << "found Global RotDataPath: "
         << RotDataProp->GlobalFilePath << "\n\n";
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

    cout << "##########################################################################################" << "\n";
    cout << "## Run - Rotation File Explorer ##########################################################" << "\n";
    cout << "##########################################################################################" << "\n\n";

    bool CheckFlag = false;

    // get global path of the RotationData file
    get_GlobalRotDataPath(Local_RotDataPath, RotDataProp);

    // count number of entries in the rotation data file
    NumberOfEntriesInRotationFile(&RotDataProp->Number_Of_Elements,
                                  RotDataProp->GlobalFilePath);

    cout << "Number of Entries: "
         << RotDataProp->Number_Of_Elements << "\n\n";

    cout << "##########################################################################################" << "\n";
    cout << "## Stop - Rotation File Explorer #########################################################" << "\n";
    cout << "##########################################################################################" << "\n\n\n\n";

    if(RotDataProp->Number_Of_Elements != 0){
        CheckFlag = true;
    }

    return CheckFlag;
}