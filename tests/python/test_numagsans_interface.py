from pathlib import Path

import pytest

from NuMagSANS import NuMagSANS


def test_generate_all_outputs_contains_expected_keys():
    outputs = NuMagSANS.generate_all_outputs()

    assert "SpinFlip_1D" in outputs
    assert "SpinFlip_2D" in outputs
    assert "Nuclear_1D" in outputs
    assert "Fourier_Gamma" in outputs
    assert len(outputs) == len(set(outputs))


def test_write_config_writes_selected_parameters(tmp_path):
    sim = NuMagSANS(executable=tmp_path / "missing_executable", workdir=tmp_path)
    config = sim.write_config(
        tmp_path / "NuMagSANSInput_test.conf",
        MagData_activate=1,
        MagDataPath="RealSpaceData/MagData",
        MagData_ReplicationImport=1,
        MagData_NumberOfReplications=8,
        RotDataLoop=1,
        RotData_User_Selection=[1, 3],
        User_Selection=[2],
        q_min=0.1,
        q_max=2.5,
        enable_outputs=["SpinFlip_1D", "Fourier_Gamma"],
    )

    text = config.read_text()
    assert "MagData_activate = 1;" in text
    assert "MagData_ReplicationImport = 1;" in text
    assert "MagData_NumberOfReplications = 8;" in text
    assert "RotDataLoop = 1;" in text
    assert "RotData_User_Selection = {1, 3};" in text
    assert "User_Selection = {2};" in text
    assert "q_min = 0.1;" in text
    assert "q_max = 2.5;" in text
    assert "SpinFlip_1D = 1;" in text
    assert "Fourier_Gamma = 1;" in text
    assert "SpinFlip_2D = 0;" in text


def test_write_config_rejects_unknown_output(tmp_path):
    sim = NuMagSANS(executable=tmp_path / "missing_executable", workdir=tmp_path)

    with pytest.raises(ValueError, match="Unknown output keys"):
        sim.write_config(tmp_path / "bad.conf", enable_outputs=["NotAnOutput"])


def test_config_clear_only_deletes_conf_by_default(tmp_path):
    sim = NuMagSANS(executable=tmp_path / "missing_executable", workdir=tmp_path)

    config = tmp_path / "input.conf"
    config.write_text("temporary config")
    sim.config_clear(config)
    assert not config.exists()

    non_config = tmp_path / "input.txt"
    non_config.write_text("not a config")
    with pytest.raises(ValueError, match="Refusing to delete non-.conf"):
        sim.config_clear(non_config)
    assert non_config.exists()

    sim.config_clear(non_config, confirm=True)
    assert not non_config.exists()


def test_run_missing_executable_raises_file_not_found(tmp_path):
    sim = NuMagSANS(executable=tmp_path / "missing_executable", workdir=tmp_path)
    config = tmp_path / "input.conf"
    config.write_text("")

    with pytest.raises(FileNotFoundError, match="Executable not found"):
        sim.run(config)


def test_constructor_accepts_explicit_paths(tmp_path):
    executable = tmp_path / "NuMagSANS"
    sim = NuMagSANS(executable=executable, workdir=tmp_path)

    assert sim.executable == executable.resolve()
    assert sim.workdir == Path(tmp_path).resolve()
