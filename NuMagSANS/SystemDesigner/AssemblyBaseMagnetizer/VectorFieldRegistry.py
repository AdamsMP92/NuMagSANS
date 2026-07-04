"""Dispatch vector-field parameter dictionaries to concrete field libraries."""

from .VectorFieldModels.OperatorKernelVectorFieldLib import unit_field as operator_kernel_unit_field
from .VectorFieldModels.SphericalVectorFieldLib import unit_field as spherical_unit_field


def _library_name(params):
    library = params.get("library", params.get("model_family"))
    if library is not None:
        return str(library).lower()

    field_type = str(params.get("field_type", "spherical")).lower()
    if field_type in {"operator_kernel", "gaussian_operator_kernel"}:
        return "operator_kernel"

    return "spherical"


def evaluate_vector_field(x, y, z, D, params):
    """Evaluate a vector field from a SystemDesigner field-parameter dictionary."""
    params = dict(params)
    library = _library_name(params)

    if library in {"spherical", "spherical_vector_field"}:
        return spherical_unit_field(x, y, z, D, params)

    if library in {"operator_kernel", "operator-kernel", "gaussian_operator_kernel"}:
        return operator_kernel_unit_field(x, y, z, D, params)

    raise ValueError(f"Unknown vector-field library '{library}'.")
