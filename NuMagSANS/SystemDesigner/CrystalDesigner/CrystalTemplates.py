import numpy as np

from .CrystalBase import crystal_base
from .CrystalOperator import crystal_centering, implicit_crystal_operator
from .CrystalRepeater import crystal_repeater


def sc_sphere_crystal(a, R, atomtype):
    """Create a spherical cut from a simple cubic one-atom lattice.

    This template is a convenience wrapper for a common construction pipeline:

        crystal base -> simple cubic repetition -> centering -> spherical cut

    It is intentionally small and explicit, so future templates can follow the
    same structure while changing only the base atoms, lattice vectors, or final
    geometric operation.

    Parameters
    ----------
    a : float
        Simple cubic lattice constant.
    R : float
        Radius of the spherical cut after centering the repeated crystal.
    atomtype : object
        Free-form atom type label for the one-atom base, for example ``"Fe"``.

    Returns
    -------
    dict
        Crystal dictionary containing the cut crystal positions, atom types,
        base mapping arrays, and lattice metadata.
    """
    # Simple cubic lattice directions.
    u1 = [1, 0, 0]
    u2 = [0, 1, 0]
    u3 = [0, 0, 1]

    # Build a one-atom crystal base. The atom type label is free-form and is not
    # chemically validated at this layer.
    base = crystal_base(
        x=0.0,
        y=0.0,
        z=0.0,
        atomtype=atomtype,
    )

    # Number of repetitions along x, y, and z before the spherical cut.
    N = 2 * np.ceil(R / a) + 1

    # Repeat the base in a simple cubic lattice.
    crystal = crystal_repeater(
        crystal_base=base,
        a1=a,
        a2=a,
        a3=a,
        N1=N,
        N2=N,
        N3=N,
        u1=u1,
        u2=u2,
        u3=u3,
    )

    # Center before cutting so the spherical shell is placed around the origin.
    crystal = crystal_centering(crystal)

    # Keep all atoms inside or on the spherical boundary:
    #     x**2 + y**2 + z**2 <= R**2
    crystal = implicit_crystal_operator(
        crystal,
        "x**2 + y**2 + z**2 - R**2",
        {"R": R},
    )

    return crystal


def sc_cylinder_crystal(a, R, h, atomtype):
    """Create a cylindrical cut from a simple cubic one-atom lattice.

    This template is a convenience wrapper for a common construction pipeline:

        crystal base -> simple cubic repetition -> centering -> cylindrical cut

    It is intentionally small and explicit, so future templates can follow the
    same structure while changing only the base atoms, lattice vectors, or final
    geometric operation.

    Parameters
    ----------
    a : float
        Simple cubic lattice constant.
    R : float
        Radius of the cylindrical cut after centering the repeated crystal.
    h : float
        height of the cylindrical cut after centering the repeated crystal.
    atomtype : object
        Free-form atom type label for the one-atom base, for example ``"Fe"``.

    Returns
    -------
    dict
        Crystal dictionary containing the cut crystal positions, atom types,
        base mapping arrays, and lattice metadata.
    """
    # Simple cubic lattice directions.
    u1 = [1, 0, 0]
    u2 = [0, 1, 0]
    u3 = [0, 0, 1]

    # Build a one-atom crystal base. The atom type label is free-form and is not
    # chemically validated at this layer.
    base = crystal_base(
        x=0.0,
        y=0.0,
        z=0.0,
        atomtype=atomtype,
    )

    # Number of repetitions along x, y, and z before the spherical cut.
    Nx = 2 * np.ceil(R / a) + 1
    Ny = 2 * np.ceil(R / a) + 1
    Nz = 2 * np.ceil(h / a) + 1

    # Repeat the base in a simple cubic lattice.
    crystal = crystal_repeater(
        crystal_base=base,
        a1=a,
        a2=a,
        a3=a,
        N1=Nx,
        N2=Ny,
        N3=Nz,
        u1=u1,
        u2=u2,
        u3=u3,
    )

    # Center before cutting so the spherical shell is placed around the origin.
    crystal = crystal_centering(crystal)

    # Keep all atoms inside or on the cylindrical boundary:
    #     x**2 + y**2 <= R**2
    crystal = implicit_crystal_operator(
        crystal,
        "x**2 + y**2 - R**2",
        {"R": R},
    )

    # Keep all atoms inside or on the cylindrical boundary:
    #    |z| <= h/2
    crystal = implicit_crystal_operator(
        crystal,
        "abs(z) - h/2",
        {"h": h},
    )

    return crystal


CRYSTAL_TEMPLATE_REGISTRY = {
    "sc_sphere_crystal": sc_sphere_crystal,
    "sc_cylinder_crystal": sc_cylinder_crystal,
}


def get_crystal_template(template_name):
    """Return a crystal template function by name."""
    if template_name not in CRYSTAL_TEMPLATE_REGISTRY:
        raise KeyError(f"Unknown crystal template: {template_name}")

    return CRYSTAL_TEMPLATE_REGISTRY[template_name]
