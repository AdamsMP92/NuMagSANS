.. toctree::
   :maxdepth: 2
   :hidden:

   Home <index>
   Getting Started <GettingStarted/index>
   Example Gallery <examples/index>
   Related Publications <RelatedPublications/index>


NuMagSANS Documentation
=======================

Overview
--------

NuMagSANS is a GPU-accelerated computational engine for the simulation
and analysis of nuclear and magnetic small-angle neutron scattering (SANS).

The aim of NuMagSANS is to provide the complete set of experimentally accessible SANS cross sections and associated observables, including unpolarized measurements as well as spin-resolved and spin-polarized channels such as SANSPOL and POLARIS. Contributions from nuclear and magnetic scattering can be treated separately or in combination, making the framework applicable beyond purely magnetic scattering problems.

NuMagSANS is designed as a deterministic single-run computational engine
with a clear execution pipeline and reproducible output.

Key Features
------------

- GPU-accelerated atomistic SANS computation
- Full polarization analysis (spin-flip, non-spin-flip, chiral terms)
- 2D and 1D scattering cross sections
- Correlation functions and pair distribution functions
- Spectral decomposition of scattering amplitudes
- Python facade for simplified workflow integration

Corresponding Publication
-------------------------

If you use NuMagSANS in your work, please cite:

**NuMagSANS: a GPU-accelerated open-source software package for the
generic computation of nuclear and magnetic small-angle neutron scattering
observables of complex systems.**

M. P. Adams and A. Michels, *Journal of Applied Crystallography* **59**,
3 (2026).

- DOI: `10.1107/S160057672600258X <https://doi.org/10.1107/S160057672600258X>`_
- Preprint: `arXiv:2601.18444 <https://doi.org/10.48550/arXiv.2601.18444>`_

.. code-block:: bibtex

   @article{adams2026numagsans_paper,
     author  = {Adams, Michael P. and Michels, Andreas},
     title   = {{\it NuMagSANS}: a GPU-accelerated open-source software package for the generic computation of nuclear and magnetic small-angle neutron scattering observables of complex systems},
     journal = {Journal of Applied Crystallography},
     year    = {2026},
     volume  = {59},
     number  = {3},
     pages   = {},
     month   = {Jun},
     doi     = {10.1107/S160057672600258X},
     url     = {https://doi.org/10.1107/S160057672600258X},
   }

License
-------

NuMagSANS is released under the MIT License.
