
.. _program_listing_file_src_OutputData_SpectralData_NuMagSANSlib_SpectralData2hdf5.h:

Program Listing for File NuMagSANSlib_SpectralData2hdf5.h
=========================================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_OutputData_SpectralData_NuMagSANSlib_SpectralData2hdf5.h>` (``src/OutputData/SpectralData/NuMagSANSlib_SpectralData2hdf5.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   
   #include "NuMagSANSlib_SpectralDataTable.h"
   
   #ifdef NUMAGSANS_ENABLE_HDF5
   #include <hdf5.h>
   #endif
   
   inline std::string build_SpectralData_hdf5_group_path(int MagData_File_Index, int StructData_File_Index = 0,
                                                         int RotData_File_Index = 0) {
       std::string group_path = "/SANS_" + std::to_string(MagData_File_Index);
   
       if (StructData_File_Index > 0) {
           group_path += "/StructData_" + std::to_string(StructData_File_Index);
       }
   
       if (RotData_File_Index > 0) {
           group_path += "/RotData_" + std::to_string(RotData_File_Index);
       }
   
       return group_path + "/AngularSpectrum";
   }
   
   #ifdef NUMAGSANS_ENABLE_HDF5
   inline void ensureSpectralHDF5Group(hid_t file_id, const std::string& group_path) {
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
   
   inline void writeSpectralHDF5Dataset(hid_t file_id, const std::string& dataset_path, unsigned long length,
                                        const float* data) {
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
   
   inline void writeSpectralHDF5Table(hid_t file_id, const std::string& table_path, unsigned long length,
                                      const std::vector<SpectralColumn>& columns) {
       ensureSpectralHDF5Group(file_id, table_path);
   
       for (const auto& column : columns) {
           writeSpectralHDF5Dataset(file_id, table_path + "/" + column.name, length, column.data);
       }
   }
   #endif
   
   inline void write2HDF5_SpectralData(InputFileData* InputData, SpectralData* SpecData, int MagData_File_Index,
                                       int StructData_File_Index = 0, int RotData_File_Index = 0) {
   #ifndef NUMAGSANS_ENABLE_HDF5
       LogSystem::write("Error: HDF5 spectral output requested, but NuMagSANS was built without HDF5 support.");
       LogSystem::write("Reconfigure with -DNUMAGSANS_ENABLE_HDF5=ON and make sure HDF5 is available.");
       throw std::runtime_error("HDF5 spectral output requested without HDF5 build support.");
   #else
       LogSystem::write("");
       LogSystem::write("write spectral decomposition data to hdf5-file...");
   
       std::filesystem::create_directories(InputData->SANSDataFoldername);
   
       std::string filename = InputData->SANSDataFoldername + "/NuMagSANS_Output.h5";
       hid_t file_id = std::filesystem::exists(filename)
                           ? H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT)
                           : H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   
       if (file_id < 0) {
           throw std::runtime_error("Could not open HDF5 output file: " + filename);
       }
   
       std::string root_path =
           build_SpectralData_hdf5_group_path(MagData_File_Index, StructData_File_Index, RotData_File_Index);
   
       unsigned int Nq = *SpecData->Nq;
       unsigned int k_max = *SpecData->k_max;
       unsigned int Nk = k_max + 1;
   
       for (const auto& comp : build_spectral_intensities(SpecData)) {
           writeSpectralHDF5Table(file_id, root_path + "/Intensities/" + comp.name, Nq,
                                  build_spectral_intensity_columns(SpecData, comp));
       }
   
       for (const auto& comp : build_spectral_amplitudes(SpecData)) {
           std::vector<float> k_values;
           writeSpectralHDF5Table(file_id, root_path + "/Amplitudes/" + comp.name, Nk,
                                  build_spectral_amplitude_columns(SpecData, comp, k_values));
       }
   
       H5Fclose(file_id);
       LogSystem::write("All spectral HDF5 datasets written successfully.");
   #endif
   }
