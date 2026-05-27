#pragma once

struct NuMagSANSData {
	ScalingFactors ScalFactors;
	NuclearData NucData, NucData_gpu;
	MagnetizationData MagData, MagData_gpu;
	StructureData StructData, StructData_gpu;
	RotationData RotData, RotData_gpu;
	ScatteringData SANSData, SANSData_gpu;
	SpectralData SpecData, SpecData_gpu;
};
