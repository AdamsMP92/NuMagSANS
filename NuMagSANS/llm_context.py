"""Version-matched usage context for AI assistants."""

from importlib.metadata import PackageNotFoundError, version
from importlib.resources import files

from .NuMagSANS import NuMagSANS


def _installed_version() -> str:
    try:
        return version("NuMagSANS")
    except PackageNotFoundError:
        return "source checkout"


def llm_context() -> str:
    """Return installed facts followed by the NuMagSANS assistant guide."""

    sim = NuMagSANS()
    guide = files("NuMagSANS").joinpath("llms-full.txt").read_text(encoding="utf-8").strip()
    conventions = ", ".join(sorted(NuMagSANS.ROT_DATA_CONVENTIONS))
    fourier_approaches = ", ".join(sorted(NuMagSANS.FOURIER_APPROACHES))
    output_formats = ", ".join(sorted(NuMagSANS.OUTPUT_FORMATS))
    outputs = ", ".join(NuMagSANS.generate_all_outputs())

    facts = "\n".join(
        [
            "# Installed NuMagSANS facts",
            "",
            f"- Python package version: {_installed_version()}",
            f"- Backend executable candidate: {sim.executable}",
            f"- Backend executable is a file: {'yes' if sim.executable.is_file() else 'no'}",
            f"- Supported Fourier approaches: {fourier_approaches}",
            f"- Accepted output format names: {output_formats}",
            "- HDF5 compiled into the backend: not detectable here; verify the CMake build setting",
            f"- Supported RotData conventions: {conventions}",
            f"- Supported output keys: {outputs}",
            "",
            "Treat these installed facts and Python introspection as authoritative for this installation.",
            "",
        ]
    )
    return facts + guide + "\n"


def main() -> None:
    """Print the usage context for terminal-based AI assistant workflows."""

    print(llm_context(), end="")
