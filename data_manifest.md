# Introduction

This file documents the  datasets included in this public repository. pCS raw data
is confidential as it's still unpublished. For this reason we generated a dummy dataset
to ensure pipeline reproducibility. pCS project names are also anonymized. 
eCS raw data is provided.


---

## Files with synthetic (dummy) data

The pCS mastersheet and SESOI file require
anonymization via synthetic data, since these are the only files containing
project-identifying information subject to confidentiality restrictions. Both
real versions (`MASTERSHEET_DECIDE_anonymized.csv`, `SESOI_DECIDE_anonymized.csv`)
are excluded from this repository; the corresponding `_DUMMY`
files below are provided.

### 1. `MASTERSHEET_DECIDE_anonymized_DUMMY.csv`

**Purpose:** Synthetic stand-in for the real `MASTERSHEET_DECIDE_anonymized.csv`,
matching its structure so the confirmatory replication assessment pipeline
(`confirmatory_replication_assessment.R`) can be run and validated.

**Rows:** One row per (project × stage × center × treatment group), following
the same exploratory/confirmatory structure as the real dataset, across 10
synthetic projects (`A`–`J`).

**Columns:**

| Column | Type | Description |
|---|---|---|
| `Project_ID` | character | Synthetic project identifier (`A`–`J`) |
| `Project_Name` | character | Same as `Project_ID` in the dummy data |
| `project_letter` | character | Same as `Project_ID` in the dummy data |
| `Stage` | character | `"exploratory"` or `"confirmatory"` |
| `Center` | character | Synthetic lab/center identifier (`Center_1`, `Center_2`, ...) |
| `Group_ES` | character | `"Ctrl"` or `"Treated"` |
| `Value` | numeric | Randomly generated raw observation (used for "raw" and "subgroup" design projects); `NA` for "precomputed" design projects |
| `Mean` | numeric | Randomly generated precomputed group mean (only populated for "precomputed" design projects) |
| `SD` | numeric | Randomly generated precomputed group SD (only populated for "precomputed" design projects) |
| `n` | integer | Precomputed group sample size (only populated for "precomputed" design projects) |
| `Subgroup` | character | Subgroup label (only populated for "subgroup" design projects, mirroring real projects with multiple measurement batches per lab) |
| `Treatment` | character | Placeholder value (`"DummyTreatment"`) for all rows |
| `planned_EU_total` | integer | Randomly generated planned total experimental units (populated on one summary row per project) |
| `planned_EU_per_group` | integer | Randomly generated planned experimental units per group (populated on one summary row per project) |

**Design across the 10 synthetic projects**: a mix of raw individual observations, precomputed
mean/SD/n only, and subgroup-level raw data requiring pooling.

**Number of confirmatory-stage laboratories per project:** varies (2–3),
matching the real dataset's range of multi-lab confirmatory designs.


### 2. `SESOI_DECIDE_anonymized_DUMMY.csv`

**Purpose:** Dataset providing the smallest-effect-size-of-interest (SESOI) value per project for the
replication assessment pipeline's "CI above SESOI" criterion.

**Rows:** One row per project.

**Columns:**

| Column | Type | Description |
|---|---|---|
| `Project_Name` | character | Synthetic project identifier (`A`–`J`), matching `Project_Name` in the mastersheet |
| `SESOI_g` | numeric | Randomly generated smallest effect size of interest (Hedges' g scale) |

---

## Other data files: provided as-is

The following files, used by `confirmatory_vs_retrospective.R` and
`protocol_comparison.R`, are included in this
repository in their real, unmodified form (no dummy/synthetic substitute).

### `mastersheet_external_studies.csv`

**Content:** Retrospective/eCS multi-laboratory replication data, sourced from
external published studies, used by `confirmatory_vs_retrospective.R`.

**Rows:** One row per (external study ID × stage × lab/center × subgroup ×
treatment group).

**Columns:**

| Column | Type | Description |
|---|---|---|
| `id` | integer | External study identifier |
| `stage` | character | `"Exploratory"` or `"Multi_lab"` |
| `center` | character | Lab/center identifier |
| `subgroup` | character | Subgroup label, where applicable (e.g., species subgroup, measurement batch) |
| `group_es` | character | `"ctrl"` or `"treated"` |
| `sample_size` | integer | Group sample size |
| `mean_difference` | numeric | Mean difference (control vs. treated) |
| `pooled_sd` | numeric | Pooled standard deviation |
| `J` | numeric | Hedges' g correction factor |
| `hedges_g` | numeric | Hedges' g effect size |
| `se_g` | numeric | Standard error of `hedges_g` |
| `EU` | integer | Experimental unit count |
| `iv_score` | numeric | Internal validity score |

### `decide_mIV_anonymized.csv`

**Content:** Per-lab minimal internal validity (mIV) scores for the pCS
projects, used by `confirmatory_vs_retrospective.R` for comparison against the
retrospective dataset.

**Rows:** One row per (project × stage) combination.

**Columns:**

| Column | Type | Description |
|---|---|---|
| `Project_ID` | character | Project/lab identifier |
| `Study_Type` | character | `"Exploratory"` or `"Confirmatory"` |
| `IV_final` | numeric | Minimal internal validity (mIV) score |

### `protocol_comparison_anonymized.csv`

**Content:** Protocol-level validity assessment scores across the four validity
domains (internal, external, statistical, translational), used by
`protocol_comparison.R` for the supervised/unsupervised heatmaps and clustering
analysis.

**Rows:** One row per (project × phase × validity domain × protocol item).

**Columns:**

| Column | Type | Description |
|---|---|---|
| `project_letter` | character | Project identifier (`A`–`K`) |
| `phase` | character | `"Exploratory Stage"` or `"Confirmatory Stage"` |
| `validity` | character | One of `"IV"`, `"EV"`, `"SV"`, `"TV"` (internal, external, statistical, translational validity) |
| `aspect` | character | Protocol sub-category within a validity domain (e.g., `"Blinding"`, `"Randomization"`) |
| `detail` | character | Specific protocol item label (e.g., `"Data analysis"`, `"Group allocation"`), mapped in analysis code to a cleaned display label (e.g., `"Blinded Analysis"`) |
| `score` | numeric | Raw item score |
| `max_score` | numeric | Maximum possible score for that item; analysis code computes `rel_score = score / max_score` |

---

## Note for users running this pipeline

`confirmatory_replication_assessment.R` reads data using the real filenames
(`MASTERSHEET_DECIDE_anonymized.csv`, `SESOI_DECIDE_anonymized.csv`), which are
not included in this repository. To run the script using the synthetic data
provided here, rename or copy the two `_DUMMY` files to match the filenames
expected by the script.

Note: Figures and tables generated this way are for pipeline-verification purposes
only and will not match the published manuscript's results, since the
underlying values in the dummy files are randomly generated.
