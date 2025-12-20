[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17933975.svg)](https://doi.org/10.5281/zenodo.17933975)

# ES-MADM III – Computational Model

This repository provides the official computational implementation of the
**ES-MADM III (Entropy Synergy Multi-Attribute Decision-Making)** framework.

The model is designed to support robust, interpretable, and numerically
well-posed multi-criteria decision analysis through entropy-based weighting,
integrated criteria importance, and stability-aware ranking mechanisms.

---

## Scope of the Repository

This repository contains **only the computational core** of the ES-MADM III model.
It is intentionally limited to numerical computations and result generation.

Visualization scripts, diagnostic plots, and presentation-level figures used in
related publications are **excluded by design**, in order to preserve clarity,
reproducibility, and methodological transparency.

---

## Repository Structure

ES-MADM-III-Computational-Model/
├── code/ # Core R implementation of ES-MADM III
├── README.md # Project description and usage notes
├── LICENSE # MIT License
---

## Case Study Data

The Excel input datasets used in the case study experiments reported in the
associated publication are provided in the `data/case_study_inputs` directory.
These datasets are directly compatible with the ES-MADM III computational core
and enable full reproducibility of the reported results.

## Requirements

- R (version ≥ 4.2)
- Standard R packages as declared in the source code

No external software or proprietary dependencies are required.

---

## Usage

The main R script located in the `code/` directory implements the complete
computational workflow of the ES-MADM III model.

Users may adapt input data and parameters according to their specific
decision-making context.

---

## Citation

If you use this software in academic work, please cite the associated Zenodo DOI
(which will be generated upon public release of this repository).

---

## License

This project is licensed under the MIT License.



