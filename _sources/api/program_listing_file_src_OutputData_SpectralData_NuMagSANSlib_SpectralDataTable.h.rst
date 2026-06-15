
.. _program_listing_file_src_OutputData_SpectralData_NuMagSANSlib_SpectralDataTable.h:

Program Listing for File NuMagSANSlib_SpectralDataTable.h
=========================================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_OutputData_SpectralData_NuMagSANSlib_SpectralDataTable.h>` (``src/OutputData/SpectralData/NuMagSANSlib_SpectralDataTable.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   
   #include <string>
   #include <vector>
   
   #include "NuMagSANSlib_SpectralDataStruct.h"
   
   struct SpectralComponent {
       std::string name;
       const float* data;
   
       SpectralComponent(const std::string& n, const float* d) : name(n), data(d) {}
   };
   
   struct SpectralColumn {
       std::string name;
       const float* data;
   
       SpectralColumn(const std::string& n, const float* d) : name(n), data(d) {}
   };
   
   inline std::vector<SpectralComponent> build_spectral_intensities(SpectralData* SpecData) {
       return {{"I_Nuc_unpolarized", SpecData->I_Nuc_unpolarized},
               {"I_Mag_unpolarized", SpecData->I_Mag_unpolarized},
               {"I_Mag_polarized", SpecData->I_Mag_polarized},
               {"I_NucMag", SpecData->I_NucMag},
               {"I_Mag_chiral", SpecData->I_Mag_chiral},
               {"I_Mag_spin_flip", SpecData->I_Mag_spin_flip},
               {"I_Mag_spin_flip_pm", SpecData->I_Mag_spin_flip_pm},
               {"I_Mag_spin_flip_mp", SpecData->I_Mag_spin_flip_mp},
               {"I_Mag_non_spin_flip_pp", SpecData->I_Mag_non_spin_flip_pp},
               {"I_Mag_non_spin_flip_mm", SpecData->I_Mag_non_spin_flip_mm},
               {"I_Mag_sanspol_p", SpecData->I_Mag_sanspol_p},
               {"I_Mag_sanspol_m", SpecData->I_Mag_sanspol_m}};
   }
   
   inline std::vector<SpectralComponent> build_spectral_amplitudes(SpectralData* SpecData) {
       return {{"A_Nuc_unpolarized", SpecData->A_Nuc_unpolarized},
               {"A_Mag_unpolarized", SpecData->A_Mag_unpolarized},
               {"A_Mag_polarized", SpecData->A_Mag_polarized},
               {"A_NucMag", SpecData->A_NucMag},
               {"A_Mag_chiral", SpecData->A_Mag_chiral},
               {"A_Mag_spin_flip", SpecData->A_Mag_spin_flip},
               {"A_Mag_spin_flip_pm", SpecData->A_Mag_spin_flip_pm},
               {"A_Mag_spin_flip_mp", SpecData->A_Mag_spin_flip_mp},
               {"A_Mag_non_spin_flip_pp", SpecData->A_Mag_non_spin_flip_pp},
               {"A_Mag_non_spin_flip_mm", SpecData->A_Mag_non_spin_flip_mm},
               {"A_Mag_sanspol_p", SpecData->A_Mag_sanspol_p},
               {"A_Mag_sanspol_m", SpecData->A_Mag_sanspol_m}};
   }
   
   inline std::vector<SpectralColumn> build_spectral_intensity_columns(SpectralData* SpecData,
                                                                       const SpectralComponent& comp) {
       unsigned int Nq = *SpecData->Nq;
       unsigned int k_max = *SpecData->k_max;
       unsigned int L = k_max + 1;
       unsigned int offset_cos = 0;
       unsigned int offset_sin = Nq * L;
   
       std::vector<SpectralColumn> columns;
       columns.emplace_back("q", SpecData->q);
   
       for (unsigned int k = 0; k <= k_max; ++k) {
           columns.emplace_back("Ic_" + std::to_string(k), comp.data + offset_cos + k * Nq);
       }
   
       for (unsigned int k = 0; k <= k_max; ++k) {
           columns.emplace_back("Is_" + std::to_string(k), comp.data + offset_sin + k * Nq);
       }
   
       return columns;
   }
   
   inline std::vector<SpectralColumn>
   build_spectral_amplitude_columns(SpectralData* SpecData, const SpectralComponent& comp, std::vector<float>& k_values) {
       unsigned int k_max = *SpecData->k_max;
       unsigned int L = k_max + 1;
   
       k_values.resize(L);
       for (unsigned int k = 0; k <= k_max; ++k) {
           k_values[k] = static_cast<float>(k);
       }
   
       std::vector<SpectralColumn> columns;
       columns.emplace_back("k", k_values.data());
       columns.emplace_back("A_cos", comp.data);
       columns.emplace_back("A_sin", comp.data + L);
   
       return columns;
   }
