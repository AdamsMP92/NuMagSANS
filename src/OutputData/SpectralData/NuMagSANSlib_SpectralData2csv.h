#pragma once

#include "NuMagSANSlib_SpectralDataTable.h"

void write_spectral_csv(const std::string& folder, const SpectralComponent& comp, SpectralData* SpecData) {
    unsigned int Nq = *SpecData->Nq;
    auto columns = build_spectral_intensity_columns(SpecData, comp);

    std::ofstream fout(folder + comp.name + ".csv");

    fout << columns.front().name;
    for (unsigned int col = 1; col < columns.size(); ++col)
        fout << "," << columns[col].name;
    fout << "\n";

    for (unsigned int i = 0; i < Nq; ++i) {
        fout << columns.front().data[i];
        for (unsigned int col = 1; col < columns.size(); ++col)
            fout << "," << columns[col].data[i];
        fout << "\n";
    }
}

void write_amplitude_csv(const std::string& folder, const SpectralComponent& comp, SpectralData* SpecData) {
    unsigned int k_max = *SpecData->k_max;
    std::vector<float> k_values;
    auto columns = build_spectral_amplitude_columns(SpecData, comp, k_values);

    std::ofstream fout(folder + comp.name + ".csv");
    if (!fout.is_open()) {
        LogSystem::write("Error opening amplitude file: " + std::string(comp.name));
        return;
    }

    fout << columns.front().name;
    for (unsigned int col = 1; col < columns.size(); ++col)
        fout << "," << columns[col].name;
    fout << "\n";

    for (unsigned int k = 0; k <= k_max; ++k) {
        fout << columns.front().data[k];
        for (unsigned int col = 1; col < columns.size(); ++col)
            fout << "," << columns[col].data[k];
        fout << "\n";
    }
}

void write2CSV_SpectralData(InputFileData* InputData, SpectralData* SpecData, int MagData_File_Index,
                            int StructData_File_Index = 0, int RotData_File_Index = 0) {
    LogSystem::write("");
    LogSystem::write("write spectral decomposition data to csv-files...");

    std::string baseFolder = InputData->SANSDataFoldername + "/SANS_" + std::to_string(MagData_File_Index) + "/";

    if (StructData_File_Index > 0) {
        baseFolder += "StructData_" + std::to_string(StructData_File_Index) + "/";
    }

    if (RotData_File_Index > 0) {
        baseFolder += "RotData_" + std::to_string(RotData_File_Index) + "/";
    }

    baseFolder += "AngularSpectrum/";

    // mkdir(baseFolder.c_str(), 0777);
    std::filesystem::create_directories(baseFolder);

    // ================================
    // Intensities
    // ================================
    std::string intensityFolder = baseFolder + "Intensities/";
    // mkdir(intensityFolder.c_str(), 0777);
    std::filesystem::create_directories(intensityFolder);

    auto intensities = build_spectral_intensities(SpecData);

    for (const auto& comp : intensities)
        write_spectral_csv(intensityFolder, comp, SpecData);

    // ================================
    // Amplitudes
    // ================================
    std::string ampFolder = baseFolder + "Amplitudes/";
    // mkdir(ampFolder.c_str(), 0777);
    std::filesystem::create_directories(ampFolder);

    auto amplitudes = build_spectral_amplitudes(SpecData);

    for (const auto& comp : amplitudes)
        write_amplitude_csv(ampFolder, comp, SpecData);

    LogSystem::write("All spectral CSV files written successfully.");
}
