from .AssemblyBase import assembly_base, constant, normal
from .AssemblyBaseWriter import write_assembly_base_objects

SC_SPHERE_CRYSTAL_TEMPLATE = "sc_sphere_crystal"


def write_monodisperse_spherical_nanoparticle_base(
    R,
    a,
    atomtype,
    n_objects,
    output_dir,
    name="Monodisperse spherical nanoparticle base",
):
    """Generate and store a monodisperse spherical nanoparticle object base.

    This is the high-level template entry point. It hides the underlying
    ``sc_sphere_crystal`` template name from examples and user scripts.
    """

    assembly_base_data = assembly_base(
        ct_temps=SC_SPHERE_CRYSTAL_TEMPLATE,
        ct_temps_params=["a", "R", "atomtype"],
        param_dist_props={
            "a": constant(a),
            "R": constant(R),
            "atomtype": constant(atomtype),
        },
        name=name,
    )

    return write_assembly_base_objects(
        assembly_base_data,
        template_name=SC_SPHERE_CRYSTAL_TEMPLATE,
        n_objects=n_objects,
        output_dir=output_dir,
        seed=None,
    )


def write_gaussian_spherical_nanoparticle_base(
    R_mean,
    R_std,
    a,
    atomtype,
    n_objects,
    output_dir,
    seed=None,
    name="Gaussian spherical nanoparticle base",
):
    """Generate and store a Gaussian spherical nanoparticle object base.

    This is the high-level template entry point. It hides the underlying
    ``sc_sphere_crystal`` template name from examples and user scripts.
    """

    assembly_base_data = assembly_base(
        ct_temps=SC_SPHERE_CRYSTAL_TEMPLATE,
        ct_temps_params=["a", "R", "atomtype"],
        param_dist_props={
            "a": constant(a),
            "R": normal(mean=R_mean, sigma=R_std),
            "atomtype": constant(atomtype),
        },
        name=name,
    )

    return write_assembly_base_objects(
        assembly_base_data,
        template_name=SC_SPHERE_CRYSTAL_TEMPLATE,
        n_objects=n_objects,
        output_dir=output_dir,
        seed=seed,
    )
