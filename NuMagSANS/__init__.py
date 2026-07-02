"""
NuMagSANS Python interface package.

Provides a high-level Python facade for executing the
NuMagSANS CUDA backend.
"""

from . import SystemDesigner
from .NuMagSANS import NuMagSANS

__all__ = ["NuMagSANS", "SystemDesigner"]
