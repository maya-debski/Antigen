# Antigen: Data Reduction Pipeline for GCMS (VIRUS-P), VIRUS-W, and VIRUS-2

[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](LICENSE)
<!--[![CI](https://github.com/grzeimann/Panacea/actions/workflows/python-tests.yml/badge.svg)](https://github.com/grzeimann/Panacea/actions/workflows/python-tests.yml)-->
<!--[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18250411.svg)](https://doi.org/10.5281/zenodo.18250411)-->

Antigen is designed to reduce data from the George and Cynthia Mitchell Spectrograph (GCMS), VIRUS-W, and the forthcoming VIRUS-2 instruments on the 2.7m Harlan J. Smith Telescope at McDonald Observatory. Antigen automates CCD reduction, wavelength and flux calibration, and the generation of science-ready datacubes and outputs FITS files and diagnostic plots at each step of the process (base CCD reduction, response measurement (optional), and cube creation). For a structured overview, see [Project documentation index](docs/index.md).

- Getting started: [Installation](docs/getting-started/installation.md) · [Quickstart](docs/getting-started/quickstart.md)
- Data products: [Overview](docs/data-products/overview.md)

## Install (conda + pip)

```bash
conda env create -f environment.yml
conda activate env_antigen
pip install . --no-deps
```

<!--Compatibility note: The environment is pinned to numpy<2.0 with astropy<6 (see [environment.yml](environment.yml) / [pyproject.toml](pyproject.toml)). If you need NumPy 2.x, upgrade Astropy to ≥6 and test locally.-->

## Quickstart (CLI)

For a step-by-step walkthrough, sample data download, and additional options, see the full Quickstart guide. It complements and expands the brief commands shown here: [docs/getting-started/quickstart.md](docs/getting-started/quickstart.md)


More examples and options: [CLI user guide](docs/user-guide/cli.md) and [Examples](docs/user-guide/examples.md).

## Data and configuration
- Required configuration files (e.g., line lists, throughput curves, config YAMLs) are bundled with the package. <!--and resolved via importlib.resources.-->
- Raw data layout: the CLI expects an specific file tree strucutre. Run `antigen_prepare_dataset.py` for Antigen compatibility.

## Troubleshooting
- Need dev commands (pre‑commit, local CI, coverage)? See [Dev quickstart](docs/development/dev-quickstart.md). For exact git workflow commands (stage/commit/push), see [Stage, commit, push](docs/development/dev-quickstart.md#stage-commit-push-exact-commands).
- Where are tests described (including tests/test_sample_data.py)? See [Tests overview](docs/development/tests-overview.md).
- More help: [FAQ](docs/faq.md).

## Citing and License
- Citation: <!--[CITATION.cff](CITATION.cff) (see also [Citation docs](docs/citation/citation.md))-->
- DOI (concept, stable): **ZENODO LINK HERE**
- DOI (version v1.0.1): **ZENODO LINK HERE**
- License: BSD‑3‑Clause ([LICENSE](LICENSE))

## Contributing
We welcome issues and pull requests. Please read [Contributing guide](docs/community/contributing.md). For a concise developer setup and CI‑equivalent checks, see [Dev quickstart](docs/development/dev-quickstart.md).

## Developers

- Maya Debski 
- Gregory Zeimann
- Jason Vestuto
