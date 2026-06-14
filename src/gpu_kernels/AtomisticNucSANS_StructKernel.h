#pragma once

// Computes purely nuclear atomistic SANS for an explicitly structured system.
// Active data layers:
//   - NuclearData: local nuclear scattering length density data.
//   - StructureData: object-center translations.
// Inactive data layers:
//   - MagnetizationData
//   - RotationData
// Use case:
//   Nuclear scattering from an assembly with explicit object positions, where
//   all objects share the same local orientation convention.
//
// Placeholder for:
// Atomistic_NucSANS_Kernel

__global__ void Atomistic_NucSANS_Kernel(NuclearData NucData, StructureData StructData, ScatteringData SANSData) {

    // Input information:
    // N     : number of atoms
    // L     : number of points in Fourier space L = N_q*N_theta
    // K     : number of particles
    // x     : x-real-space position data in units of nano-meters
    // y     : y-real-space position data in units of nano-meters
    // z     : z-real-space position data in units of nano-meters
    // qy    : qy-Fourier-space coordinate in units of inverse nano-meters
    // qz    : qz-Fourier-space coordinate in units of inverse nano-meters
    // theta : theta-angle in Fourier space [theta = arctan2(qz, qy)] in radiant

    // output information:

    unsigned long int L = (*SANSData.N_q) * (*SANSData.N_theta);
    unsigned long int N_avg = *NucData.N_avg;
    unsigned long int K = *NucData.K;
    unsigned long int W = *NucData.TotalAtomNumber;

    // float v = 1.0/((float)  (*NucData.K)) * powf(1.0/((float) (*NucData.N)), 2); // pre factor
    float v = 1.0 / ((float)W) * 1.0 / ((float)N_avg);

    int i = blockIdx.x * blockDim.x + threadIdx.x;

    float nuc_real = 0.0;
    float nuc_imag = 0.0;

    float Nuc_real = 0.0;
    float Nuc_imag = 0.0;

    float Y = 0.0;
    float Z = 0.0;

    float Psi = 0.0;

    float cos_val = 0.0;
    float sin_val = 0.0;

    unsigned long int N_cum = 0;
    unsigned long int N_act = 0;

    if (i < L) {
        for (int k = 0; k < K; k++) {

            nuc_real = 0.0;
            nuc_imag = 0.0;

            N_cum = NucData.N_cum[k];
            N_act = NucData.N_act[k];

            for (int l = 0; l < N_act; l++) {
                // atomic position composition
                // X = MagData.RotMat[0] * (MagData.x[l+k*N] + StructData.x[k]) \
                //  + MagData.RotMat[3] * (MagData.y[l+k*N] + StructData.y[k]) \
                // + MagData.RotMat[6] * (MagData.z[l+k*N] + StructData.z[k]);
                Y = NucData.RotMat[1] * (NucData.x[l + N_cum] + StructData.x[k]) +
                    NucData.RotMat[4] * (NucData.y[l + N_cum] + StructData.y[k]) +
                    NucData.RotMat[7] * (NucData.z[l + N_cum] + StructData.z[k]);
                Z = NucData.RotMat[2] * (NucData.x[l + N_cum] + StructData.x[k]) +
                    NucData.RotMat[5] * (NucData.y[l + N_cum] + StructData.y[k]) +
                    NucData.RotMat[8] * (NucData.z[l + N_cum] + StructData.z[k]);

                // phase function
                Psi = Y * SANSData.qy_2D[i] + Z * SANSData.qz_2D[i];

                // cosine and sine values
                cos_val = cosf(Psi);
                sin_val = sinf(Psi);

                nuc_real += NucData.Nuc[l + N_cum] * cos_val;
                nuc_imag -= NucData.Nuc[l + N_cum] * sin_val;
            }

            Nuc_real += nuc_real;
            Nuc_imag += nuc_imag;
        }

        // nuclear SANS cross section projected in (qz, qy)-plane
        SANSData.S_Nuc_2D_unpolarized[i] = v * (Nuc_real * Nuc_real + Nuc_imag * Nuc_imag);
    }
}
