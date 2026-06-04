// File         : NuMagSANSlib_Directory.h
// Author       : Michael Philipp ADAMS, M.Sc.
// Company      : University of Luxembourg
// Department   : Department of Physics and Materials Sciences
// Group        : NanoMagnetism Group
// Group Leader : Prof. Andreas Michels
// Version      : 22 November 2024
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


int count_NumberOfFolders(string Directory){

     struct dirent* entity;
     DIR* dir = opendir(Directory.c_str());           // open directory stored in tmp_char
     entity = readdir(dir);              // readdir command reads file and foldernames in directory
     string subfolder_name;
     int counter = 0;         // counter that counts how many folders are in the MagData folder
     while(entity != NULL){
              subfolder_name = entity->d_name;
         if(subfolder_name != "." && subfolder_name != ".." && subfolder_name != ".DS_Store" && entity->d_type == DT_DIR){
             counter += 1;
         }
         entity = readdir(dir);
     }
     closedir(dir);
     return counter;

 }

int count_NumberOfElements(string Directory){

     struct dirent* entity;
     DIR* dir = opendir(Directory.c_str());           // open directory stored in tmp_char
     entity = readdir(dir);              // readdir command reads file and foldernames in directory
     string element_name;
     int counter = 0;         // counter that counts how many folders are in the MagData folder
     while(entity != NULL){
         element_name = entity->d_name;
         if(element_name != "." && element_name != ".." && element_name != ".DS_Store"){
             counter += 1;
         }
         entity = readdir(dir);
     }
     closedir(dir);
     return counter;

 }


 void read_FolderNames(string Directory, string*FolderNames, int Number_Of_FolderNames){

     struct dirent* entity;
     DIR* dir = opendir(Directory.c_str());           // open directory stored in tmp_char
     entity = readdir(dir);              // readdir command reads file and foldernames in directory
     string subfolder_name;
     int counter = 0;         // counter that counts how many folders are in the MagData folder
     while(entity != NULL && counter < Number_Of_FolderNames){
         subfolder_name = entity->d_name;
         if(subfolder_name != "." && subfolder_name != ".." && subfolder_name != ".DS_Store"){
             FolderNames[counter] = subfolder_name;
             counter += 1;
         }
         entity = readdir(dir);
     }
     closedir(dir);

}



