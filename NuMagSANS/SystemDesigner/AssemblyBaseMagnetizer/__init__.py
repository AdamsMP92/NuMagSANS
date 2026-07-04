"""Magnetization tools for materialized assembly-base objects."""

from .AssemblyBaseMagnetizer import write_spherical_magnetization_data
from .MagnetizationBase import magnetization_base, sample_field_params
from .MagnetizationBaseTemplates import (
    spherical_random_chirality_vortex_magnetization_base,
    spherical_vortex_magnetization_base,
)
from .VectorFieldModels import (
    GaussianKernel,
    RealSpaceComponentSpec,
    SechKernel,
    alpha_profile,
    build_operator_kernel_model,
    spherical_unit_field,
)
from .VectorFieldRegistry import evaluate_vector_field

unit_field = spherical_unit_field

__all__ = [
    "GaussianKernel",
    "RealSpaceComponentSpec",
    "SechKernel",
    "alpha_profile",
    "build_operator_kernel_model",
    "evaluate_vector_field",
    "magnetization_base",
    "sample_field_params",
    "spherical_random_chirality_vortex_magnetization_base",
    "spherical_vortex_magnetization_base",
    "unit_field",
    "write_spherical_magnetization_data",
]
