# DECIDE: Replication Assessment of Preclinical Multi-Laboratory Studies

This repository contains the analysis code for the DECIDE meta-research project, which evaluates the replicability of preclinical multi-laboratory studies and benchmarks statistical criteria for assessing replication success.

## Project Overview

The DECIDE project (Decision-Enabling Confirmation of Innovative Discoveries and Exploratory Evidence) accompanied twelve preclinical confirmatory multi-laboratory studies funded by the German Federal Ministry of Research, Technology and Space (BMFTR). This repository provides:

- **Confirmatory replication assessment**: effect size estimation, meta-analysis, and replication criteria applied to DECIDE confirmatory studies
- **Retrospective comparison**: comparison with published preclinical multi-laboratory studies from the literature
- **Protocol comparison**: analysis of internal, external, translational, and statistical validity across exploratory and confirmatory stages
- **Simulation framework**: stress-testing of replication criteria under realistic preclinical conditions

## Repository Structure

```
DECIDE/
├── Main.R                                    # Master script — runs full pipeline
├── packages.R                                 # All package dependencies
├── all_functions.R                            # Shared utility functions
├── confirmatory_replication_assessment.R      # DECIDE confirmatory analysis
├── confirmatory_vs_retrospective.R            # Comparison with retrospective dataset
├── protocol_comparison.R                      # Validity scoring and radar plots
├── simulation_replication_criteria.R          # Simulation framework
├── data/                                      # Input data (see Data Availability)
├── results/                                   # Output figures and tables
├── simulation_results/                        # Simulation outputs per dataset
├── renv.lock                                  # Package version lockfile
└── README.md
```

## Reproducibility

This project uses [`renv`](https://rstudio.github.io/renv/) to manage package versions.

### Setup

1. Clone the repository
2. Open `DECIDE.Rproj` in RStudio — **always open via the `.Rproj` file**
3. Install exact package versions:
```r
renv::restore()
```
4. Place data files in the `/data` folder (see Data Availability below)
5. Run the full analysis pipeline:
```r
source("Main.R")
```

> **Note:** Scripts must be run in order via `Main.R`. Running individual scripts out of order may fail due to dependencies between scripts.

> **Note on data:** Because the primary confirmatory study (pCS/DECIDE) data are unpublished, this repository includes a **synthetic (randomly generated) replacement** for the pCS mastersheet and SESOI files, matching their structure exactly (`MASTERSHEET_DECIDE_anonymized_DUMMY.csv`, `SESOI_DECIDE_anonymized_DUMMY.csv`; see `data_manifest.md`). This allows the confirmatory replication assessment pipeline to be run end-to-end and verified for reproducibility, but figures and results generated from this synthetic data are for pipeline-verification purposes only and will **not** match those reported in the manuscript.

### R Version
R 4.5.1 — see `renv.lock` for exact package versions.

## Data Availability

| File | Description |
|------|-------------|
| `MASTERSHEET_DECIDE_anonymized_DUMMY.csv` | Synthetic replacement for the pCS/DECIDE raw animal study mastersheet, matching its structure (see `data_manifest.md`) |
| `SESOI_DECIDE_anonymized_DUMMY.csv` | Synthetic replacement for pCS smallest effect sizes of interest per project (see `data_manifest.md`) |
| `decide_mIV_anonymized.csv` | Minimal internal validity scores for DECIDE (pCS) projects |
| `mastersheet_external_studies.csv` | Retrospective (eCS) multi-laboratory dataset |
| `protocol_comparison_anonymized.csv` | Protocol validity scoring data (internal, external, statistical, translational) |
| `empyrical_effect_size_datasets/` | Empirical effect size distributions for the simulation framework |

**Note on pCS raw data:** The primary confirmatory study (pCS/DECIDE) raw dataset is **not publicly available**, as the underlying studies are unpublished. Aggregated/summary data sufficient to verify the manuscript's reported results are available from the corresponding author upon reasonable request.

Place all files in the `/data` folder before running the analysis.

## Simulation Framework

The simulation stress-tests replication criteria under realistic preclinical conditions by modelling a two-stage research trajectory:

- **True effects** sampled from empirical distributions (Bonapersona et al. 2021; Carneiro et al. 2018; Rosso et al. 2022)
- **Exploratory studies**: n = 5, 10, 15, or 20 per group; only significant results (p < 0.05) advance to confirmatory stage
- **Confirmatory studies**: k = 3 laboratories, each powered at 80% based on the observed exploratory effect (max n = 50 per group), pooled via fixed-effects meta-analysis
- **Shrinkage**: confirmatory true effect = exploratory true effect × (1 − s), where s ∈ {0, 0.2, 0.5, 0.8, 0.99}

Performance is evaluated using:
- **Success rate** across effect size classes, shrinkage levels, and exploratory sample sizes
- **False positive rate** at s = 0.99 (near-null confirmatory effect)
- **Precision and recall** with low shrinkage (s ∈ {0, 0.2}) as ground truth positives and s = 0.99 as ground truth negatives
- **Shrinkage sensitivity**: correlation between shrinkage level and criterion success rate

## Citation

> Pellegrini P. et al. (2026). *Title TBA*. Journal TBA. DOI: TBA

Code DOI: [Zenodo — to be added]

## License

MIT License — see `LICENSE.md` for details.

## Contact

Pasquale Pellegrini\
Berlin Institute of Health at Charité (BIH-QUEST)\
pasquale.pellegrini@bih-charite.de
