
.. _program_listing_file_src_NuMagSANSlib_ReadWrite.h:

Program Listing for File NuMagSANSlib_ReadWrite.h
=================================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_NuMagSANSlib_ReadWrite.h>` (``src/NuMagSANSlib_ReadWrite.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   // File         : NuMagSANSlib_ReadWrite.h
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
   
   using namespace std;
   
   // Read and Write Function library
   
   // ##############################################################################################################################################
   // ##############################################################################################################################################
   
   void WriteCSV_double_vector(unsigned long int N, double *data, string filename){
       ofstream fout;
       fout.open(filename);
       for(unsigned long int n=0; n<N; n++) fout << data[n] << "\n";
       fout.close();
   }
   
   void WriteCSV_float_vector(unsigned long int N, float *data, string filename){
       ofstream fout;
       fout.open(filename);
       for(unsigned long int n=0; n<N; n++) fout << data[n] << "\n";
       fout.close();
   }
   
   void ReadCSV_double_vector(unsigned long int N, double *data, string filename){
       ifstream fin;
       fin.open(filename);
       for(unsigned long int n=0; n<N; n++) fin >> data[n];
       fin.close();
   }
   
   void ReadCSV_float_vector(unsigned long int N, float *data, string filename){
       ifstream fin;
       fin.open(filename);
       for(unsigned long int n=0; n<N; n++) fin >> data[n];
       fin.close();
   }
   
   void ReadCSV_double_matrix(unsigned long int N, unsigned long int M, double **data, string filename){
       ifstream fin;
       fin.open(filename);
       for(unsigned long int n=0; n<N; n++){
           for(unsigned long int m=0; m<M; m++) fin >> data[n][m];
       }
       fin.close();
   }
   
   void WriteCSV_float_matrix(unsigned long int N, unsigned long int M, float **data, string filename){
       ofstream fout;
       fout.open(filename);
       for(unsigned long int n=0; n<N; n++){
           for(unsigned long int m=0; m<M; m++){
               fout << data[n][m];
               if(m<M-1) fout << ",";
           }
           fout << "\n";
       }
       fout.close();
   }
   
   
   
   
   
   
   
   
   /*
   
   void NumberOfEntriesInNucFile(int *NumberOfColumns, string filename){
   
       ifstream fin;
       fin.open(filename);
       std::string line;
   
       float x, y, z, Nuc;
   
       int line_counter = 0;
       int error_counter = 0;
   
       while(std::getline(fin, line)){
           std::istringstream ss(line);
           if(ss >> x >> y >> z >> Nuc){
               line_counter += 1;          
           } else{
               error_counter ++;
   //          std::cerr << "Error in row: " << line_counter + error_counter << ": " << line << "\n"; 
           }
       }
       fin.close();
       *NumberOfColumns = line_counter;
   
   }
   */
   
   /*
   void NumberOfNonZeroMagneticMomentsInFile(int *NumberOfNonZeroMoments, int *NumberOfColumns, string filename){
   
       ifstream fin;
       fin.open(filename);
       std::string line;
   
       float x, y, z, mx, my, mz;
   
       int line_counter = 0;
       int moment_counter = 0;
       int error_counter = 0;
   
       while(std::getline(fin, line)){
           std::istringstream ss(line);
           if(ss >> x >> y >> z >> mx >> my >> mz){
               if(mx != 0.0 || my != 0.0 || mz != 0.0){
                   moment_counter += 1;
               }
               line_counter += 1;          
           } else{
               error_counter ++;
               std::cerr << "Error in row: " << line_counter + error_counter << ": " << line << "\n"; 
           }
       }
       fin.close();
       *NumberOfColumns = line_counter;
       *NumberOfNonZeroMoments = moment_counter;
   
   }
   */
   
   /*
   void CountColumnsAndRowsInMagDataFile(int *Number_Of_Rows, int *Number_Of_Columns, string filename){
   
       ifstream fin;
       fin.open(filename);
       string line;
       //float ghost_buf = 0.0;
   
       *Number_Of_Rows = 0;
       *Number_Of_Columns = 0;
       while(getline(fin, line)){
           //cout << line.size() << "\n";
           *Number_Of_Rows += 1;
       }
   }
   */
