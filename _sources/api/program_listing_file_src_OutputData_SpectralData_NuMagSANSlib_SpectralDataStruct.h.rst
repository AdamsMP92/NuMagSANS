
.. _program_listing_file_src_OutputData_SpectralData_NuMagSANSlib_SpectralDataStruct.h:

Program Listing for File NuMagSANSlib_SpectralDataStruct.h
==========================================================

|exhale_lsh| :ref:`Return to documentation for file <file_src_OutputData_SpectralData_NuMagSANSlib_SpectralDataStruct.h>` (``src/OutputData/SpectralData/NuMagSANSlib_SpectralDataStruct.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   
   struct SpectralData {
       unsigned int* Nq;     // number of q-values
       unsigned int* Ntheta; // number of theta-values
       unsigned int* k_max;  // maximum number of modes
       float* dtheta;        // angular step-size
       float* q;
   
       // the arrays are single column and contain the angular sine and cosine spectra
       // the total array length is 2 * Nq * (k_max + 1)
       float* I_Nuc_unpolarized;      // nuclear SANS cross section
       float* I_Mag_unpolarized;      // unpolarized magnetic SANS cross section
       float* I_Mag_polarized;        // polarized magnetic SANS cross section
       float* I_NucMag;               // nuclear-magnetic interference SANS cross section
       float* I_Mag_spin_flip;        // spin-flip magnetic SANS cross section
       float* I_Mag_chiral;           // chiral magnetic SANS cross section
       float* I_Mag_spin_flip_pm;     // pm-spin-flip magnetic SANS cross section
       float* I_Mag_spin_flip_mp;     // mp-spin-flip magnetic SANS cross section
       float* I_Mag_non_spin_flip_pp; // pp-non-spin-flip magnetic SANS cross section
       float* I_Mag_non_spin_flip_mm; // mm-non-spin-flip magnetic SANS cross section
       float* I_Mag_sanspol_p;        // p-sanspol magnetic SANS cross section
       float* I_Mag_sanspol_m;        // m-sanspol magnetic SANS cross section
   
       // normalized angular amplitude spectrum
       // the total array length is 2 * (k_max + 1)
       float* A_Nuc_unpolarized;      // nuclear SANS cross section
       float* A_Mag_unpolarized;      // unpolarized magnetic SANS cross section
       float* A_Mag_polarized;        // polarized magnetic SANS cross section
       float* A_NucMag;               // nuclear-magnetic interference SANS cross section
       float* A_Mag_spin_flip;        // spin-flip magnetic SANS cross section
       float* A_Mag_chiral;           // chiral magnetic SANS cross section
       float* A_Mag_spin_flip_pm;     // pm-spin-flip magnetic SANS cross section
       float* A_Mag_spin_flip_mp;     // mp-spin-flip magnetic SANS cross section
       float* A_Mag_non_spin_flip_pp; // pp-non-spin-flip magnetic SANS cross section
       float* A_Mag_non_spin_flip_mm; // mm-non-spin-flip magnetic SANS cross section
       float* A_Mag_sanspol_p;        // p-sanspol magnetic SANS cross section
       float* A_Mag_sanspol_m;        // m-sanspol magnetic SANS cross section
   };
