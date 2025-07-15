# INSTALL

Instructions on how to install the `antigen` package

## Summary Install

- (optional) `./install_conda.sh` 
- `conda env create -f environment.yml`
- `conda activate env_antigen`
- `pip install . --no-deps`

## Conda Overview 

- `conda` is strongly recommended for environment and dependency management with this package.
- In examples here, all external package dependencies are managed by `conda`
- The reduction package itself is as a `pip` package which is installed into the `conda` env with `pip install`

## Using Conda

- This repo is organized as a pip package which can be installed in many different types of environments.
- Antigen installation examples here recommend and assume use of a conda environment.
- For scripted installation of conda on macOS and Linux, see [`./install_conda.sh`](`./install_conda.sh`)
- For manual installation of conda
    - Download conda installer from https://github.com/conda-forge/miniforgeInstall 
    - Install conda under `$HOME/opt/conda` for user not system 
    - Review [`./install_conda.sh`](`./install_conda.sh`) for initial conda configuration steps

## Installing into the Conda Environment

- After conda is installed ...
- Create the conda environment from the `environemnt.yml`:
    - Create a conda environment for antigen: `conda env create -f environment.yml`
    - Activate conda env: `conda activate env_antigen`
- Pip Install the antigen package into conda env: 
    - `cd <repo-root-path>`
    - `pip install . --no-deps`
    - Optional install methods:
        - for testing during development: `pip install --editable . --no-deps`
        - for testing dependencies after install: `pip install --no-build-isolation --no-deps .`

## Details: Pip Packaging, for Development and Deployment

Building distribution packages:
- activate the conda env: `conda activate env_antigen`
- install the pip package `build` tool: `conda install -c conda-forge build`
- build a src dist: 
    - build: `python -m build --sdist`
    - verify: `ls -l ./dist/antigen-2025.07.15.tar.gz`
    - install: `pip install ./dist/antigen-2025.07.15.tar.gz --no-deps`
- build a wheel: 
    - build: `python -m build --wheel`
    - verify: `ls -l ./dist/antigen-2025.07.15-py3-none-any.whl`
    - install: `pip install ./dist/antigen-2025.07.15-py3-none-any.whl --no-deps`
