Tutorial Slides
===============

The tutorial slide decks are published as PDF files. They are built from the
slide sources in ``docs/tutorials`` and copied to the root ``tutorials`` folder
with ``make publish``.

Each PDF title page contains a version date generated at build time.

The tutorial scripts generate their data with the core installation. Install the
tutorial dependencies when you also want the diagnostic vector-field plots:

.. code-block:: bash

   pip install -e ".[tutorials]"

.. list-table::
   :header-rows: 1
   :widths: 20 50 30

   * - Tutorial
     - Topic
     - PDF
   * - Tutorial 1
     - HPC installation
     - `Open PDF <../tutorial_1/slides.pdf>`__
   * - Tutorial 2
     - Getting started workflow
     - `Open PDF <../tutorial_2/slides.pdf>`__
   * - Tutorial 3
     - Vector-field workflows
     - `Open PDF <../tutorial_3/slides.pdf>`__
   * - Tutorial 4
     - Helical magnet workflow
     - `Open PDF <../tutorial_4/slides.pdf>`__
   * - Tutorial 5
     - Longitudinal helix sweep
     - `Open PDF <../tutorial_5/slides.pdf>`__

Build Notes
-----------

The published PDFs are served by Sphinx through ``html_extra_path`` in
``docs/conf.py``. This copies the root ``tutorials`` folder into the generated
HTML site during the documentation build.

When a slide deck changes, rebuild and publish it from its source folder:

.. code-block:: bash

   make slides.pdf
   make publish

or use:

.. code-block:: bash

   make final
