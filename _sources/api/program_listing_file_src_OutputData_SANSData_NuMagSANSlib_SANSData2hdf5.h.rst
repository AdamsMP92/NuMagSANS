
.. _program_listing_file_src_OutputData_SANSData_NuMagSANSlib_SANSData2hdf5.h:

Program Listing for File NuMagSANSlib_SANSData2hdf5.h
=====================================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_OutputData_SANSData_NuMagSANSlib_SANSData2hdf5.h>` (``src/OutputData/SANSData/NuMagSANSlib_SANSData2hdf5.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   
   #include "NuMagSANSlib_SANSDataTable.h"
   
   #ifdef NUMAGSANS_ENABLE_HDF5
   #include <hdf5.h>
   #endif
   
   inline std::string build_SANSData_hdf5_group_path(int MagData_File_Index, int StructData_File_Index = 0,
                                                     int RotData_File_Index = 0) {
       std::string group_path = "/SANS_" + std::to_string(MagData_File_Index);
   
       if (StructData_File_Index > 0) {
           group_path += "/StructData_" + std::to_string(StructData_File_Index);
       }
   
       if (RotData_File_Index > 0) {
           group_path += "/RotData_" + std::to_string(RotData_File_Index);
       }
   
       return group_path;
   }
   
   #ifdef NUMAGSANS_ENABLE_HDF5
   inline void ensureHDF5Group(hid_t file_id, const std::string& group_path) {
       std::string current_path;
       std::stringstream ss(group_path);
       std::string token;
   
       while (std::getline(ss, token, '/')) {
           if (token.empty()) {
               continue;
           }
   
           current_path += "/" + token;
   
           if (H5Lexists(file_id, current_path.c_str(), H5P_DEFAULT) <= 0) {
               hid_t group_id = H5Gcreate2(file_id, current_path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
               if (group_id < 0) {
                   throw std::runtime_error("Could not create HDF5 group: " + current_path);
               }
               H5Gclose(group_id);
           }
       }
   }
   
   inline void writeHDF5Dataset(hid_t file_id, const std::string& dataset_path, unsigned long length, const float* data) {
       hsize_t dims[1] = {static_cast<hsize_t>(length)};
       hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
       if (dataspace_id < 0) {
           throw std::runtime_error("Could not create HDF5 dataspace for: " + dataset_path);
       }
   
       if (H5Lexists(file_id, dataset_path.c_str(), H5P_DEFAULT) > 0) {
           H5Ldelete(file_id, dataset_path.c_str(), H5P_DEFAULT);
       }
   
       hid_t dataset_id = H5Dcreate2(file_id, dataset_path.c_str(), H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT,
                                     H5P_DEFAULT, H5P_DEFAULT);
       if (dataset_id < 0) {
           H5Sclose(dataspace_id);
           throw std::runtime_error("Could not create HDF5 dataset: " + dataset_path);
       }
   
       if (H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0) {
           H5Dclose(dataset_id);
           H5Sclose(dataspace_id);
           throw std::runtime_error("Could not write HDF5 dataset: " + dataset_path);
       }
   
       H5Dclose(dataset_id);
       H5Sclose(dataspace_id);
   }
   
   inline void writeHDF5Table(hid_t file_id, const std::string& table_path, unsigned long length,
                              const std::vector<Column>& columns) {
       ensureHDF5Group(file_id, table_path);
   
       for (const auto& column : columns) {
           writeHDF5Dataset(file_id, table_path + "/" + column.name, length, column.data);
       }
   }
   #endif
   
   inline void write2HDF5table_ScatteringData(InputFileData* InputData, ScatteringData* SANSData, int MagData_File_Index,
                                              int StructData_File_Index = 0, int RotData_File_Index = 0) {
   #ifndef NUMAGSANS_ENABLE_HDF5
       LogSystem::write("Error: HDF5 output requested, but NuMagSANS was built without HDF5 support.");
       LogSystem::write("Reconfigure with -DNUMAGSANS_ENABLE_HDF5=ON and make sure HDF5 is available.");
       throw std::runtime_error("HDF5 output requested without HDF5 build support.");
   #else
       LogSystem::write("");
       LogSystem::write("write scattering data to hdf5-file...");
   
       std::filesystem::create_directories(InputData->SANSDataFoldername);
   
       std::string filename = InputData->SANSDataFoldername + "/NuMagSANS_Output.h5";
       hid_t file_id = std::filesystem::exists(filename)
                           ? H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT)
                           : H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   
       if (file_id < 0) {
           throw std::runtime_error("Could not open HDF5 output file: " + filename);
       }
   
       std::string root_path =
           build_SANSData_hdf5_group_path(MagData_File_Index, StructData_File_Index, RotData_File_Index);
   
       if (InputData->output_fourier_correlation_matrix_flag || any_active(InputData->OutFlags.SANS2D)) {
           unsigned long L = (*SANSData->N_q) * (*SANSData->N_theta);
           writeHDF5Table(file_id, root_path + "/SANS2D", L, build_SANS2D_columns(InputData, SANSData));
           LogSystem::write("SANS2D hdf5 output finished...");
       }
   
       if (any_active(InputData->OutFlags.SANS1D)) {
           writeHDF5Table(file_id, root_path + "/SANS1D", (*SANSData->N_q), build_SANS1D_columns(InputData, SANSData));
           LogSystem::write("SANS1D hdf5 output finished...");
       }
   
       if (any_active(InputData->OutFlags.Corr1D) || any_active(InputData->OutFlags.PairDist1D)) {
           writeHDF5Table(file_id, root_path + "/Corr1D", (*SANSData->N_r), build_Corr1D_columns(InputData, SANSData));
           LogSystem::write("Corr1D hdf5 output finished...");
       }
   
       if (any_active(InputData->OutFlags.Corr2D)) {
           unsigned long L = (*SANSData->N_r) * (*SANSData->N_alpha);
           writeHDF5Table(file_id, root_path + "/Corr2D", L, build_Corr2D_columns(InputData, SANSData));
           LogSystem::write("Corr2D hdf5 output finished...");
       }
   
       H5Fclose(file_id);
   #endif
   }
