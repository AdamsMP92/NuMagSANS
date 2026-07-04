"""High-level templates for local magnetization parameter bases."""

from ..AssemblyBaseDesigner.AssemblyBase import constant, uniform
from .MagnetizationBase import magnetization_base

SPHERICAL_VECTOR_FIELD_TEMPLATE = "spherical_unit_field"


def spherical_vortex_magnetization_base(
    kappa=1.0,
    profile_type="linear",
    xi_type="cylindrical_xi",
    turns=1.0,
    chirality=1.0,
    name="Spherical vortex magnetization base",
):
    """Create a local spherical-vortex magnetization parameter base.

    The returned base only describes local field parameters. It does not sample
    global object rotations or translations.
    """
    return magnetization_base(
        field_templates=SPHERICAL_VECTOR_FIELD_TEMPLATE,
        field_template_params=[
            "field_type",
            "profile_type",
            "xi_type",
            "kappa",
            "turns",
            "chirality",
        ],
        param_dist_props={
            "field_type": constant("vortex"),
            "profile_type": constant(profile_type),
            "xi_type": constant(xi_type),
            "kappa": constant(kappa),
            "turns": constant(turns),
            "chirality": constant(chirality),
        },
        name=name,
    )


def spherical_random_chirality_vortex_magnetization_base(
    kappa=1.0,
    profile_type="linear",
    xi_type="cylindrical_xi",
    turns=1.0,
    name="Spherical random-chirality vortex magnetization base",
):
    """Create a local vortex base with chirality sampled uniformly in [-1, 1)."""
    return magnetization_base(
        field_templates=SPHERICAL_VECTOR_FIELD_TEMPLATE,
        field_template_params=[
            "field_type",
            "profile_type",
            "xi_type",
            "kappa",
            "turns",
            "chirality",
        ],
        param_dist_props={
            "field_type": constant("vortex"),
            "profile_type": constant(profile_type),
            "xi_type": constant(xi_type),
            "kappa": constant(kappa),
            "turns": constant(turns),
            "chirality": uniform(-1.0, 1.0),
        },
        name=name,
    )
