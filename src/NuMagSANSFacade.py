from pathlib import Path
from typing import Iterable
import subprocess


class NuMagSANSFacade:

    # ------------------------------------------------------------
    # Output naming structure (mirrors C++ logic)
    # ------------------------------------------------------------

    PREFIXES = [
        "Nuclear",
        "Unpolarized",
        "Polarized",
        "NuclearMagnetic",
        "SpinFlip",
        "Chiral",
        "PM_SpinFlip",
        "MP_SpinFlip",
        "PP_NonSpinFlip",
        "MM_NonSpinFlip",
        "P_SANSPOL",
        "M_SANSPOL",
    ]

    SUFFIXES = [
        "_2D",
        "_1D",
        "_Corr_2D",
        "_Corr_1D",
        "_PairDist_1D",
    ]

    SPECIAL_OUTPUTS = ["Fourier_Gamma"]

    @classmethod
    def generate_all_outputs(cls):
        outputs = []
        for prefix in cls.PREFIXES:
            for suffix in cls.SUFFIXES:
                outputs.append(prefix + suffix)
        outputs.extend(cls.SPECIAL_OUTPUTS)
        return outputs

    # ------------------------------------------------------------

    def __init__(
        self,
        executable: str | Path = "./NuMagSANS",
        workdir: str | Path | None = None,
    ):
        self.executable = Path(executable).resolve()
        self.workdir = Path(workdir).resolve() if workdir else Path.cwd()

    # ------------------------------------------------------------
    # CONFIG GENERATION
    # ------------------------------------------------------------

    def write_config(
        self,
        filename: str | Path,
        *,
        # foldernames
        NucDataPath="RealSpaceData/NucData",
        MagDataPath="RealSpaceData/MagData",
        StructDataFilename="RealSpaceData/StructData.csv",
        foldernameSANSData="NuMagSANS_Output",

        # selection
        NucData_activate=0,
        MagData_activate=0,
        StructData_activate=0,
        Exclude_Zero_Moments=0,

        # fourier
        Fourier_Approach="atomistic",

        # loop
        Loop_Modus=0,
        Loop_From=1,
        Loop_To=20,
        User_Selection=[1],

        # units
        XYZ_Unit_Factor=1,

        # micromagnetics
        Cell_Nuclear_SLD=8e14,
        Cell_Magnetization=486e3,
        Cuboid_Cell_Size=(2, 2, 2),

        # scattering
        Scattering_Volume_V=2.618e-24,

        # rotation
        RotMat_alpha=0.0,
        RotMat_beta=0.0,

        # polarization
        Polarization=(0.0, 0.0, 1.0),

        # q / r
        Number_Of_q_Points=1000,
        Number_Of_theta_Points=1000,
        Number_Of_r_Points=1000,
        Number_Of_alpha_Points=1000,
        q_max=3.0,
        r_max=15.0,

        # angular spectrum
        k_max=10,
        Angular_Spec=0,

        # outputs
        enable_outputs: Iterable[str] = (),
    ):

        filename = Path(filename)
        enabled = set(enable_outputs)
        all_outputs = set(self.generate_all_outputs())

        unknown = enabled - all_outputs
        if unknown:
            raise ValueError(f"Unknown output keys: {unknown}")

        with open(filename, "w") as f:

            def W(key, val):
                f.write(f"{key} = {val};\n")

            # ---------------------------
            # foldernames
            # ---------------------------
            W("NucDataPath", NucDataPath)
            W("MagDataPath", MagDataPath)
            W("StructDataFilename", StructDataFilename)
            W("foldernameSANSData", foldernameSANSData)

            # ---------------------------
            # selection
            # ---------------------------
            W("NucData_activate", NucData_activate)
            W("MagData_activate", MagData_activate)
            W("StructData_activate", StructData_activate)
            W("Exclude_Zero_Moments", Exclude_Zero_Moments)

            # ---------------------------
            # fourier
            # ---------------------------
            W("Fourier_Approach", Fourier_Approach)

            # ---------------------------
            # loop
            # ---------------------------
            W("Loop_Modus", Loop_Modus)
            W("Loop_From", Loop_From)
            W("Loop_To", Loop_To)
            W("User_Selection", "{" + ", ".join(map(str, User_Selection)) + "}")

            # ---------------------------
            # units
            # ---------------------------
            W("XYZ_Unit_Factor", XYZ_Unit_Factor)

            # ---------------------------
            # micromagnetics
            # ---------------------------
            W("Cell_Nuclear_SLD", Cell_Nuclear_SLD)
            W("Cell_Magnetization", Cell_Magnetization)
            W("Cuboid_Cell_Size_x", Cuboid_Cell_Size[0])
            W("Cuboid_Cell_Size_y", Cuboid_Cell_Size[1])
            W("Cuboid_Cell_Size_z", Cuboid_Cell_Size[2])

            # ---------------------------
            # scattering
            # ---------------------------
            W("Scattering_Volume_V", Scattering_Volume_V)

            # ---------------------------
            # rotation
            # ---------------------------
            W("RotMat_alpha", RotMat_alpha)
            W("RotMat_beta", RotMat_beta)

            # ---------------------------
            # polarization
            # ---------------------------
            W("Polarization_x", Polarization[0])
            W("Polarization_y", Polarization[1])
            W("Polarization_z", Polarization[2])

            # ---------------------------
            # q / r
            # ---------------------------
            W("Number_Of_q_Points", Number_Of_q_Points)
            W("Number_Of_theta_Points", Number_Of_theta_Points)
            W("Number_Of_r_Points", Number_Of_r_Points)
            W("Number_Of_alpha_Points", Number_Of_alpha_Points)
            W("q_max", q_max)
            W("r_max", r_max)

            # ---------------------------
            # angular spectrum
            # ---------------------------
            W("k_max", k_max)
            W("Angular_Spec", Angular_Spec)

            # ---------------------------
            # outputs
            # ---------------------------
            for key in sorted(all_outputs):
                W(key, 1 if key in enabled else 0)

        return filename


    # ------------------------------------------------------------
    # Clear Config File
    # ------------------------------------------------------------

    def config_clear(self, config: str | Path, *, confirm: bool = False):
        """
        Delete a configuration file.
        If confirm=False, only deletes files ending in '.conf'.
        """
        config = Path(config)

        if not config.exists():
            return

        if not confirm and config.suffix != ".conf":
            raise ValueError("Refusing to delete non-.conf file without confirm=True")

        config.unlink()

    # ------------------------------------------------------------
    # EXECUTION
    # ------------------------------------------------------------

    def run(self, config: str | Path, *, check: bool = True):

        if not self.executable.exists():
            raise FileNotFoundError(
                f"Executable not found: {self.executable}"
            )

        config = Path(config).resolve()

        return subprocess.run(
            [str(self.executable), str(config)],
            cwd=self.workdir,
            check=check,
        )
