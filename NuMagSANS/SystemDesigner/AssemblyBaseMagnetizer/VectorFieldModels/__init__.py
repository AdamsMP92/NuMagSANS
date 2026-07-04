"""Vector-field model families used by the SystemDesigner magnetizer."""

from .OperatorKernelVectorFieldLib import (
    GaussianKernel,
    RealSpaceComponentSpec,
    SechKernel,
    build_operator_kernel_model,
)
from .OperatorKernelVectorFieldLib import (
    unit_field as operator_kernel_unit_field,
)
from .SphericalVectorFieldLib import alpha_profile
from .SphericalVectorFieldLib import unit_field as spherical_unit_field

__all__ = [
    "GaussianKernel",
    "RealSpaceComponentSpec",
    "SechKernel",
    "alpha_profile",
    "build_operator_kernel_model",
    "operator_kernel_unit_field",
    "spherical_unit_field",
]
