#pragma once

inline void ExportDataMulti(InputFileData* InputData, ScalingFactors* ScalFactors, ScatteringData* SANSData,
                            ScatteringData* SANSData_gpu, SpectralData* SpecData, SpectralData* SpecData_gpu,
                            int Data_File_Index, int StructData_File_Index = 0, int RotData_File_Index = 0) {

    // copy scattering data from GPU to RAM ###################################################
    copyGPU2RAM_ScatteringData(SANSData, SANSData_gpu);

    // scaling of the scattering data on RAM ##################################################
    LogSystem::write("scaling of SANSdata...");
    scale_ScatteringData(ScalFactors, SANSData, InputData);

    // write scattering data ##################################################################
    if (InputData->SANSData_Output_Format == "csv") {
        write2CSVtable_ScatteringData(InputData, SANSData, Data_File_Index, StructData_File_Index, RotData_File_Index);
    } else if (InputData->SANSData_Output_Format == "hdf5") {
        write2HDF5table_ScatteringData(InputData, SANSData, Data_File_Index, StructData_File_Index, RotData_File_Index);
    }

    if (InputData->AngularSpec_activate_flag) {
        // copy spectral data from GPU to RAM #################################################
        copyGPU2RAM_SpectralData(SpecData, SpecData_gpu);

        // scaling of the spectral data on RAM ################################################
        scale_SpectralData(ScalFactors, SpecData, InputData);

        // write spectral data to csv files ###################################################
        write2CSV_SpectralData(InputData, SpecData, Data_File_Index, StructData_File_Index, RotData_File_Index);
    }
}

inline void ExportData(InputFileData* InputData, NuMagSANSData* Data, int Data_File_Index,
                       int StructData_File_Index = 0, int RotData_File_Index = 0) {

    ExportDataMulti(InputData, &Data->ScalFactors, &Data->SANSData, &Data->SANSData_gpu, &Data->SpecData,
                    &Data->SpecData_gpu, Data_File_Index, StructData_File_Index, RotData_File_Index);
}
