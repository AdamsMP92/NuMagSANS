from pathlib import Path

import numpy as np
import pytest

from NuMagSANS.SystemDesigner.Workflows import write_spherical_replication_vectorfield_sweep


def test_spherical_replication_workflow_writes_magnetic_and_rotation_data(tmp_path):
    summary = write_spherical_replication_vectorfield_sweep(
        output_dir=tmp_path,
        R=3.0,
        a=1.0,
        atomtype="Fe",
        n_replications=8,
        field_parameter_cases=[
            {
                "field_type": "vortex",
                "profile_type": "linear",
                "xi_type": "cylindrical_xi",
                "kappa": 1.0,
            },
            {
                "library": "operator_kernel",
                "kernel_type": "gaussian",
                "sigma": [0.35, 0.35, 0.35],
                "component_specs": ["-a * dy", "a * dx", "c"],
                "operator_parameters": {"a": 0.8, "c": 0.55},
                "normalize": True,
            },
        ],
        n_rotdata=2,
        rotation_seed=3,
    )

    real_space_dir = Path(summary["real_space_dir"])
    field_files = [Path(filename) for filename in summary["field_summary"]["output_files"]]
    rot_files = [Path(filename) for filename in summary["rotdata_summary"]["output_files"]]

    assert real_space_dir == tmp_path / "RealSpaceData"
    assert len(field_files) == 2
    assert len(rot_files) == 2
    assert all(filename.exists() for filename in field_files)
    assert all(filename.exists() for filename in rot_files)

    for filename in field_files:
        data = np.loadtxt(filename)
        assert data.ndim == 2
        assert data.shape[1] == 6

    for filename in rot_files:
        data = np.loadtxt(filename)
        assert data.shape == (8, 3)

    config_hints = summary["config_hints"]
    assert config_hints["MagData_ReplicationImport"] == 1
    assert config_hints["MagData_NumberOfReplications"] == 8
    assert config_hints["RotDataLoop"] == 1
    assert config_hints["RotDataLoop_To"] == 2
    assert config_hints["Loop_Modus"] == 1
    assert config_hints["Loop_To"] == 2


def test_spherical_replication_workflow_validates_inputs(tmp_path):
    with pytest.raises(ValueError, match="n_replications"):
        write_spherical_replication_vectorfield_sweep(
            output_dir=tmp_path,
            R=3.0,
            a=1.0,
            atomtype="Fe",
            n_replications=0,
        )

    with pytest.raises(ValueError, match="n_rotdata"):
        write_spherical_replication_vectorfield_sweep(
            output_dir=tmp_path,
            R=3.0,
            a=1.0,
            atomtype="Fe",
            n_replications=1,
            n_rotdata=0,
        )

    with pytest.raises(ValueError, match="beta_min"):
        write_spherical_replication_vectorfield_sweep(
            output_dir=tmp_path,
            R=3.0,
            a=1.0,
            atomtype="Fe",
            n_replications=1,
            n_rotdata=2,
            beta_min=[0.0],
        )

    with pytest.raises(KeyError, match="field_type"):
        write_spherical_replication_vectorfield_sweep(
            output_dir=tmp_path,
            R=3.0,
            a=1.0,
            atomtype="Fe",
            n_replications=1,
            field_parameter_cases=[{"profile_type": "linear"}],
        )
