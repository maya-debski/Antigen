# Installation

This page covers local installation and dependencies.

## Dependencies
- Python 3.9+
- ipython, NumPy, SciPy, Astropy, Matplotlib, PyYAML, scikit-learn, Pandas, psutil, seaborn, photutils, getcalspec, astroplan, astroscrappy

## Create a conda environment and install
```bash
conda env create -f environment.yml  # run from the repository root
conda activate env_antigen
# Install Antigen (use --no-deps because conda provided most packages)
pip install . --no-deps
```
See also: [Quickstart](../getting-started/quickstart.md)