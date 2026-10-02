import subprocess
import sys
from pathlib import Path

from NuMagSANS import llm_context

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def test_llm_context_contains_live_facts_and_canonical_guide():
    context = llm_context()

    assert context.startswith("# Installed NuMagSANS facts")
    assert "Supported RotData conventions:" in context
    assert "xyz" in context
    assert "Supported output keys:" in context
    assert "SpinFlip_1D" in context
    assert "# Using NuMagSANS — guide for AI assistants" in context
    assert "The incoming neutron beam is along the Cartesian x-axis" in context


def test_python_module_prints_llm_context():
    result = subprocess.run(
        [sys.executable, "-m", "NuMagSANS"],
        cwd=REPOSITORY_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    assert result.stderr == ""
    assert result.stdout.startswith("# Installed NuMagSANS facts")
    assert "Backend executable candidate:" in result.stdout


def test_llms_index_links_published_context_files():
    index = (REPOSITORY_ROOT / "llms.txt").read_text(encoding="utf-8")

    assert "https://adamsmp92.github.io/NuMagSANS/llms-full.txt" in index
    assert "AIAssistedSimulations.html" in index
