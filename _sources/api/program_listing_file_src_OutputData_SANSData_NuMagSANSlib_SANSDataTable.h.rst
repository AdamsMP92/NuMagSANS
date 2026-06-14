
.. _program_listing_file_src_OutputData_SANSData_NuMagSANSlib_SANSDataTable.h:

Program Listing for File NuMagSANSlib_SANSDataTable.h
=====================================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_OutputData_SANSData_NuMagSANSlib_SANSDataTable.h>` (``src/OutputData/SANSData/NuMagSANSlib_SANSDataTable.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   
   #include <string>
   #include <vector>
   
   #include "NuMagSANSlib_SANSDataStruct.h"
   
   struct Column {
       std::string name;
       const float* data;
       Column(const std::string& n, const float* d) : name(n), data(d) {}
   };
   
   std::vector<Column> build_SANS2D_columns(InputFileData* InputData, ScatteringData* SANSData) {
       std::vector<Column> columns;
   
       columns.emplace_back("qz", SANSData->qz_2D);
       columns.emplace_back("qy", SANSData->qy_2D);
       columns.emplace_back("q", SANSData->q_2D);
       columns.emplace_back("theta", SANSData->theta_2D);
   
       if (InputData->output_fourier_correlation_matrix_flag) {
           columns.emplace_back("Gxx_real", SANSData->Gxx_real);
           columns.emplace_back("Gxx_imag", SANSData->Gxx_imag);
           columns.emplace_back("Gyy_real", SANSData->Gyy_real);
           columns.emplace_back("Gyy_imag", SANSData->Gyy_imag);
           columns.emplace_back("Gzz_real", SANSData->Gzz_real);
           columns.emplace_back("Gzz_imag", SANSData->Gzz_imag);
           columns.emplace_back("Gxy_real", SANSData->Gxy_real);
           columns.emplace_back("Gxy_imag", SANSData->Gxy_imag);
           columns.emplace_back("Gyx_real", SANSData->Gyx_real);
           columns.emplace_back("Gyx_imag", SANSData->Gyx_imag);
           columns.emplace_back("Gxz_real", SANSData->Gxz_real);
           columns.emplace_back("Gxz_imag", SANSData->Gxz_imag);
           columns.emplace_back("Gzx_real", SANSData->Gzx_real);
           columns.emplace_back("Gzx_imag", SANSData->Gzx_imag);
           columns.emplace_back("Gyz_real", SANSData->Gyz_real);
           columns.emplace_back("Gyz_imag", SANSData->Gyz_imag);
           columns.emplace_back("Gzy_real", SANSData->Gzy_real);
           columns.emplace_back("Gzy_imag", SANSData->Gzy_imag);
       }
   
       auto& f = InputData->OutFlags.SANS2D;
   
       if (f.Nuclear)
           columns.emplace_back("S_N", SANSData->S_Nuc_2D_unpolarized);
       if (f.Unpolarized)
           columns.emplace_back("S_M", SANSData->S_Mag_2D_unpolarized);
       if (f.NuclearMagnetic)
           columns.emplace_back("S_NM", SANSData->S_NucMag_2D);
       if (f.Polarized)
           columns.emplace_back("S_P", SANSData->S_Mag_2D_polarized);
       if (f.Chiral)
           columns.emplace_back("S_chi", SANSData->S_Mag_2D_chiral);
       if (f.SpinFlip)
           columns.emplace_back("S_sf", SANSData->S_Mag_2D_spin_flip);
       if (f.PM_SpinFlip)
           columns.emplace_back("S_pm", SANSData->S_Mag_2D_spin_flip_pm);
       if (f.MP_SpinFlip)
           columns.emplace_back("S_mp", SANSData->S_Mag_2D_spin_flip_mp);
       if (f.PP_NonSpinFlip)
           columns.emplace_back("S_pp", SANSData->S_Mag_2D_non_spin_flip_pp);
       if (f.MM_NonSpinFlip)
           columns.emplace_back("S_mm", SANSData->S_Mag_2D_non_spin_flip_mm);
       if (f.P_SANSPOL)
           columns.emplace_back("S_p", SANSData->S_Mag_2D_sanspol_p);
       if (f.M_SANSPOL)
           columns.emplace_back("S_m", SANSData->S_Mag_2D_sanspol_m);
   
       return columns;
   }
   
   std::vector<Column> build_SANS1D_columns(InputFileData* InputData, ScatteringData* SANSData) {
       std::vector<Column> columns;
   
       columns.emplace_back("q", SANSData->q_1D);
   
       auto& f = InputData->OutFlags.SANS1D;
   
       if (f.Nuclear)
           columns.emplace_back("I_N", SANSData->S_Nuc_1D_unpolarized);
       if (f.Unpolarized)
           columns.emplace_back("I_M", SANSData->S_Mag_1D_unpolarized);
       if (f.NuclearMagnetic)
           columns.emplace_back("I_NM", SANSData->S_NucMag_1D);
       if (f.Polarized)
           columns.emplace_back("I_P", SANSData->S_Mag_1D_polarized);
       if (f.Chiral)
           columns.emplace_back("I_chi", SANSData->S_Mag_1D_chiral);
       if (f.SpinFlip)
           columns.emplace_back("I_sf", SANSData->S_Mag_1D_spin_flip);
       if (f.PM_SpinFlip)
           columns.emplace_back("I_pm", SANSData->S_Mag_1D_spin_flip_pm);
       if (f.MP_SpinFlip)
           columns.emplace_back("I_mp", SANSData->S_Mag_1D_spin_flip_mp);
       if (f.PP_NonSpinFlip)
           columns.emplace_back("I_pp", SANSData->S_Mag_1D_non_spin_flip_pp);
       if (f.MM_NonSpinFlip)
           columns.emplace_back("I_mm", SANSData->S_Mag_1D_non_spin_flip_mm);
       if (f.P_SANSPOL)
           columns.emplace_back("I_p", SANSData->S_Mag_1D_sanspol_p);
       if (f.M_SANSPOL)
           columns.emplace_back("I_m", SANSData->S_Mag_1D_sanspol_m);
   
       return columns;
   }
   
   std::vector<Column> build_Corr1D_columns(InputFileData* InputData, ScatteringData* SANSData) {
       std::vector<Column> columns;
   
       columns.emplace_back("r", SANSData->r_1D);
   
       auto& p = InputData->OutFlags.PairDist1D;
       auto& c = InputData->OutFlags.Corr1D;
   
       if (p.Nuclear)
           columns.emplace_back("p_N", SANSData->p_Nuc_unpolarized);
       if (p.Unpolarized)
           columns.emplace_back("p_M", SANSData->p_Mag_unpolarized);
       if (p.NuclearMagnetic)
           columns.emplace_back("p_NM", SANSData->p_NucMag);
       if (p.Polarized)
           columns.emplace_back("p_P", SANSData->p_Mag_polarized);
       if (p.Chiral)
           columns.emplace_back("p_chi", SANSData->p_Mag_chiral);
       if (p.SpinFlip)
           columns.emplace_back("p_sf", SANSData->p_Mag_spin_flip);
       if (p.PM_SpinFlip)
           columns.emplace_back("p_pm", SANSData->p_Mag_spin_flip_pm);
       if (p.MP_SpinFlip)
           columns.emplace_back("p_mp", SANSData->p_Mag_spin_flip_mp);
       if (p.PP_NonSpinFlip)
           columns.emplace_back("p_pp", SANSData->p_Mag_non_spin_flip_pp);
       if (p.MM_NonSpinFlip)
           columns.emplace_back("p_mm", SANSData->p_Mag_non_spin_flip_mm);
       if (p.P_SANSPOL)
           columns.emplace_back("p_p", SANSData->p_Mag_sanspol_p);
       if (p.M_SANSPOL)
           columns.emplace_back("p_m", SANSData->p_Mag_sanspol_m);
   
       if (c.Nuclear)
           columns.emplace_back("c_N", SANSData->c_Nuc_unpolarized);
       if (c.Unpolarized)
           columns.emplace_back("c_M", SANSData->c_Mag_unpolarized);
       if (c.NuclearMagnetic)
           columns.emplace_back("c_NM", SANSData->c_NucMag);
       if (c.Polarized)
           columns.emplace_back("c_P", SANSData->c_Mag_polarized);
       if (c.Chiral)
           columns.emplace_back("c_chi", SANSData->c_Mag_chiral);
       if (c.SpinFlip)
           columns.emplace_back("c_sf", SANSData->c_Mag_spin_flip);
       if (c.PM_SpinFlip)
           columns.emplace_back("c_pm", SANSData->c_Mag_spin_flip_pm);
       if (c.MP_SpinFlip)
           columns.emplace_back("c_mp", SANSData->c_Mag_spin_flip_mp);
       if (c.PP_NonSpinFlip)
           columns.emplace_back("c_pp", SANSData->c_Mag_non_spin_flip_pp);
       if (c.MM_NonSpinFlip)
           columns.emplace_back("c_mm", SANSData->c_Mag_non_spin_flip_mm);
       if (c.P_SANSPOL)
           columns.emplace_back("c_p", SANSData->c_Mag_sanspol_p);
       if (c.M_SANSPOL)
           columns.emplace_back("c_m", SANSData->c_Mag_sanspol_m);
   
       return columns;
   }
   
   std::vector<Column> build_Corr2D_columns(InputFileData* InputData, ScatteringData* SANSData) {
       std::vector<Column> columns;
   
       columns.emplace_back("rz", SANSData->rz_2D);
       columns.emplace_back("ry", SANSData->ry_2D);
       columns.emplace_back("r", SANSData->r_2D);
       columns.emplace_back("alpha", SANSData->alpha_2D);
   
       auto& f = InputData->OutFlags.Corr2D;
   
       if (f.Nuclear)
           columns.emplace_back("C_N", SANSData->Corr_Nuc_2D_unpolarized);
       if (f.Unpolarized)
           columns.emplace_back("C_M", SANSData->Corr_Mag_2D_unpolarized);
       if (f.NuclearMagnetic)
           columns.emplace_back("C_NM", SANSData->Corr_NucMag_2D);
       if (f.Polarized)
           columns.emplace_back("C_P", SANSData->Corr_Mag_2D_polarized);
       if (f.Chiral)
           columns.emplace_back("C_chi", SANSData->Corr_Mag_2D_chiral);
       if (f.SpinFlip)
           columns.emplace_back("C_sf", SANSData->Corr_Mag_2D_spin_flip);
       if (f.PM_SpinFlip)
           columns.emplace_back("C_pm", SANSData->Corr_Mag_2D_spin_flip_pm);
       if (f.MP_SpinFlip)
           columns.emplace_back("C_mp", SANSData->Corr_Mag_2D_spin_flip_mp);
       if (f.PP_NonSpinFlip)
           columns.emplace_back("C_pp", SANSData->Corr_Mag_2D_non_spin_flip_pp);
       if (f.MM_NonSpinFlip)
           columns.emplace_back("C_mm", SANSData->Corr_Mag_2D_non_spin_flip_mm);
       if (f.P_SANSPOL)
           columns.emplace_back("C_p", SANSData->Corr_Mag_2D_sanspol_p);
       if (f.M_SANSPOL)
           columns.emplace_back("C_m", SANSData->Corr_Mag_2D_sanspol_m);
   
       return columns;
   }
