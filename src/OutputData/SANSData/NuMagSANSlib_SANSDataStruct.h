#pragma once

struct ScatteringData {

    unsigned int* N_q;     // number of q-values
    unsigned int* N_theta; // number of theta-values
    unsigned int* N_r;     // number of r-values
    unsigned int* N_alpha; // number of alpha-values

    float* Polarization; // polarization vector [Px, Py, Pz]

    float* q_min; // minimum q-value
    float* q_max; // maximum q-value
    float* r_max; // maximum r-value

    float* dq;     // step-size of q
    float* dtheta; // step-size of theta
    float* dr;     // step-size of r
    float* dalpha; // step-size of alpha

    float* qy_2D;    // 2D qy-scattering vector
    float* qz_2D;    // 2D qz-scattering vector
    float* q_2D;     // 2D q-scattering vector
    float* theta_2D; // 2D theta-scattering vector

    float* ry_2D;    // 2D ry-correlation space coordinate
    float* rz_2D;    // 2D rz-correlation space coordinate
    float* r_2D;     // 2D r-correlation space coordinate
    float* alpha_2D; // 2D angle-correlation space coordinate

    float* q_1D; // 1D scattering vector
    float* r_1D; // 1D correlation space coordinate

    float* Gxx_real; // real-part xx-magnetic Fourier correlation function
    float* Gyy_real; // real-part yy-magnetic Fourier correlation function
    float* Gzz_real; // real part zz-magnetic Fourier correlation function
    float* Gxy_real; // real part xy-magnetic Fourier correlation function
    float* Gyx_real; // real part yx-magnetic Fourier correlation function
    float* Gxz_real; // real part xz-magnetic Fourier correlation function
    float* Gzx_real; // real part zx-magnetic Fourier correlation function
    float* Gyz_real; // real part yz-magnetic Fourier correlation function
    float* Gzy_real; // real part zy-magnetic Fourier correlation function

    float* Gxx_imag; // real-part xx-magnetic Fourier correlation function
    float* Gyy_imag; // real-part yy-magnetic Fourier correlation function
    float* Gzz_imag; // real-part zz-magnetic Fourier correlation function
    float* Gxy_imag; // real-part xy-magnetic Fourier correlation function
    float* Gyx_imag; // real-part yx-magnetic Fourier correlation function
    float* Gxz_imag; // real-part xz-magnetic Fourier correlation function
    float* Gzx_imag; // real-part zx-magnetic Fourier correlation function
    float* Gyz_imag; // real-part yz-magnetic Fourier correlation function
    float* Gzy_imag; // real-part zy-magnetic Fourier correlation function

    float* S_Nuc_2D_unpolarized;      // nuclear SANS cross section
    float* S_Mag_2D_unpolarized;      // unpolarized magnetic SANS cross section
    float* S_Mag_2D_polarized;        // polarized magnetic SANS cross section
    float* S_NucMag_2D;               // nuclear-magnetic interference SANS cross section
    float* S_Mag_2D_spin_flip;        // spin-flip magnetic SANS cross section
    float* S_Mag_2D_chiral;           // chiral magnetic SANS cross section
    float* S_Mag_2D_spin_flip_pm;     // pm-spin-flip magnetic SANS cross section
    float* S_Mag_2D_spin_flip_mp;     // mp-spin-flip magnetic SANS cross section
    float* S_Mag_2D_non_spin_flip_pp; // pp-non-spin-flip magnetic SANS cross section
    float* S_Mag_2D_non_spin_flip_mm; // mm-non-spin-flip magnetic SANS cross section
    float* S_Mag_2D_sanspol_p;        // p-sanspol magnetic SANS cross section
    float* S_Mag_2D_sanspol_m;        // m-sanspol magnetic SANS cross section

    float* Corr_Nuc_2D_unpolarized;
    float* Corr_Mag_2D_unpolarized;
    float* Corr_Mag_2D_polarized;
    float* Corr_NucMag_2D;
    float* Corr_Mag_2D_spin_flip;
    float* Corr_Mag_2D_chiral;
    float* Corr_Mag_2D_spin_flip_pm;
    float* Corr_Mag_2D_spin_flip_mp;
    float* Corr_Mag_2D_non_spin_flip_pp;
    float* Corr_Mag_2D_non_spin_flip_mm;
    float* Corr_Mag_2D_sanspol_p;
    float* Corr_Mag_2D_sanspol_m;

    float* S_Nuc_1D_unpolarized;
    float* S_Mag_1D_unpolarized;
    float* S_Mag_1D_polarized;
    float* S_NucMag_1D;
    float* S_Mag_1D_chiral;
    float* S_Mag_1D_spin_flip;
    float* S_Mag_1D_spin_flip_pm;
    float* S_Mag_1D_spin_flip_mp;
    float* S_Mag_1D_non_spin_flip_pp;
    float* S_Mag_1D_non_spin_flip_mm;
    float* S_Mag_1D_sanspol_p;
    float* S_Mag_1D_sanspol_m;

    float* p_Nuc_unpolarized;
    float* p_Mag_unpolarized;
    float* p_Mag_polarized;
    float* p_NucMag;
    float* p_Mag_chiral;
    float* p_Mag_spin_flip;
    float* p_Mag_spin_flip_pm;
    float* p_Mag_spin_flip_mp;
    float* p_Mag_non_spin_flip_pp;
    float* p_Mag_non_spin_flip_mm;
    float* p_Mag_sanspol_p;
    float* p_Mag_sanspol_m;

    float* c_Nuc_unpolarized;
    float* c_Mag_unpolarized;
    float* c_Mag_polarized;
    float* c_NucMag;
    float* c_Mag_chiral;
    float* c_Mag_spin_flip;
    float* c_Mag_spin_flip_pm;
    float* c_Mag_spin_flip_mp;
    float* c_Mag_non_spin_flip_pp;
    float* c_Mag_non_spin_flip_mm;
    float* c_Mag_sanspol_p;
    float* c_Mag_sanspol_m;
};
