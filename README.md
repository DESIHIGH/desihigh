# DESI High: School of the Dark Universe

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/DESIHIGH/desihigh/HEAD?urlpath=lab/tree/notebooks)
[![GitHub](https://img.shields.io/github/license/DESIHIGH/desihigh)](LICENSE.rst)

Welcome to DESI High ! This repository contains the code, ressources, notebooks and data used in the DESI High program, an initiative to teach high school students about the Dark Energy Spectrocopic Instrument ([DESI]((https://www.desi.lbl.gov/))) and cosmology in general.

🚨 For more information about the DESI High program, or access to the notebooks, please visit the [DESI High website](https://desihigh.github.io/desihigh/).

## Installation
To install the `desihigh` package, you can use pip using the following command:

```bash
pip install git+https://github.com/DESIHIGH/desihigh.git
```

You can also clone the repository and install it locally:

```bash
git clone https://github.com/DESIHIGH/desihigh.git
cd desihigh
pip install .
```

### Usage
The `desihigh` package provides various modules and functions used in the notebooks and exercises of the DESI High program.

## Repository Structure
The repository is structured as follows:

```
desihigh/
├── desihigh/                # Main package directory
├── notebooks/               # Jupyter notebooks for exercises and lessons
|   ├── English/             # Notebooks in English
|   │   ├── cheatsheets/     # Cheat sheets for quick reference
|   │   └── <notebooks>
|   ├── .../                 # Notebooks in other languages (e.g., French, Spanish)
├── data/                    # Sample data files used in the notebooks
├── images/                  # Images used in the documentation and notebooks
├── attic/                   # Archive of old notebooks and resources
├── environment.yml          # Python package dependencies for Binder
├── LICENSE.rst              # License file
├── pyproject.toml           # Project metadata and dependencies
└── README.md                # This file
```

## Contributing
We welcome contributions to the DESI High project! If you have suggestions, improvements, or bug fixes, please open an issue or submit a pull request.

### Contact
You can get in touch with the project contributors via the [DESI High website](https://desihigh.github.io/desihigh/) or via the [forum](https://github.com/DESIHIGH/desihigh/discussions)