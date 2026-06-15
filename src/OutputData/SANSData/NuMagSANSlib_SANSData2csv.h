#pragma once

#include "NuMagSANSlib_SANSDataTable.h"

void writeCSV(const std::string& filename, unsigned long length, const std::vector<Column>& columns) {
    std::ofstream fout(filename);
    if (!fout.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    // Header
    for (size_t i = 0; i < columns.size(); ++i) {
        fout << columns[i].name;
        if (i < columns.size() - 1)
            fout << ",";
    }
    fout << "\n";

    // Data
    for (unsigned long n = 0; n < length; ++n) {
        for (size_t i = 0; i < columns.size(); ++i) {
            fout << columns[i].data[n];
            if (i < columns.size() - 1)
                fout << ",";
        }
        fout << "\n";
    }

    fout.close();
}

void write2CSVtable_ScatteringData(InputFileData* InputData, ScatteringData* SANSData, int MagData_File_Index,
                                   int StructData_File_Index = 0, int RotData_File_Index = 0) {
    LogSystem::write("");
    LogSystem::write("write scattering data to csv-files...");

    std::string target_foldername = InputData->SANSDataFoldername + "/SANS_" + std::to_string(MagData_File_Index) + "/";

    if (StructData_File_Index > 0) {
        target_foldername += "StructData_" + std::to_string(StructData_File_Index) + "/";
    }

    if (RotData_File_Index > 0) {
        target_foldername += "RotData_" + std::to_string(RotData_File_Index) + "/";
    }

    // mkdir(target_foldername.c_str(), 0777);
    std::filesystem::create_directories(target_foldername);

    // ------------------- SANS2D -------------------
    if (InputData->output_fourier_correlation_matrix_flag || any_active(InputData->OutFlags.SANS2D)) {
        unsigned long L = (*SANSData->N_q) * (*SANSData->N_theta);

        auto columns = build_SANS2D_columns(InputData, SANSData);

        writeCSV(target_foldername + "SANS2D.csv", L, columns);

        LogSystem::write("SANS2D finished...");
    }

    // ------------------- SANS1D -------------------
    if (any_active(InputData->OutFlags.SANS1D)) {
        auto columns = build_SANS1D_columns(InputData, SANSData);

        writeCSV(target_foldername + "SANS1D.csv", (*SANSData->N_q), columns);

        LogSystem::write("SANS1D finished...");
    }

    // ------------------- Corr1D -------------------
    if (any_active(InputData->OutFlags.Corr1D) || any_active(InputData->OutFlags.PairDist1D)) {
        auto columns = build_Corr1D_columns(InputData, SANSData);

        writeCSV(target_foldername + "Corr1D.csv", (*SANSData->N_r), columns);

        LogSystem::write("Corr1D finished...");
    }

    // ------------------- Corr2D -------------------
    if (any_active(InputData->OutFlags.Corr2D)) {
        unsigned long L = (*SANSData->N_r) * (*SANSData->N_alpha);

        auto columns = build_Corr2D_columns(InputData, SANSData);

        writeCSV(target_foldername + "Corr2D.csv", L, columns);

        LogSystem::write("Corr2D finished...");
    }
}
