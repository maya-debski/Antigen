#!/bin/bash
# Install conda from conda-forge download.

# Detect platform and architecture
OS="$(uname -s)"
ARCH="$(uname -m)"

# Normalize architecture
case "$ARCH" in
    x86_64)
        ARCH="x86_64"
        ;;
    arm64|aarch64)
        ARCH="aarch64"
        ;;
    *)
        echo "Unsupported architecture: $ARCH"
        exit 1
        ;;
esac

# Normalize OS
case "$OS" in
    Linux)
        PLATFORM="Linux"
        ;;
    Darwin)
        PLATFORM="MacOSX"
        ;;
    *)
        echo "Unsupported OS: $OS"
        exit 1
        ;;
esac

# Construct the installer filename
INSTALLER="Miniforge3-${PLATFORM}-${ARCH}.sh"
URL="https://github.com/conda-forge/miniforge/releases/latest/download/${INSTALLER}"

# Download the installer
echo "Downloading ${INSTALLER} from ${URL}..."
curl -LO "$URL"

# Make executable
chmod +x "$INSTALLER"

# Run installer in silent/batch mode
echo "[CONDA] Running conda installer..."
bash $INSTALLER -b -f -p $HOME/opt/conda

# Initialize the conda shell config:
echo "[CONDA] Initializing conda install..."
source $HOME/opt/conda/etc/profile.d/conda.sh
conda init bash  # update user ~/.bashrc so you don't have to run etc/profile.d/conda.sh every time
conda activate base

# Conda config
conda config --add channels conda-forge
conda config --add channels defaults
conda config --set auto_activate_base false  # Disable activation of conda `base` env

# Update packages
echo "[CONDA] Updating conda packages..."
conda activate base
conda update -n base -c conda-forge conda
conda update --all
conda clean --all

# test activate
echo "[CONDA] test activating, switch to conda python env"
conda activate base   # to activate base conda env
which python          # should show conda install path

# test deactivate
echo "[CONDA] test deactivate, switch to system python env"
conda deactivate  # to deactivate conda, switch back to system paths
which python      # should show system python path

# Done!
echo "[CONDA] Install complete."
