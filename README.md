# Pollen limitation, flower visitors and human disturbance in *Buddleja delavayi*

Data and R code for a two-year field study (2025–2026) of the threatened shrub *Buddleja delavayi* (Scrophulariaceae) in Jianchuan County, Yunnan, China. The study combines a pollen-supplementation experiment, flower-visitor observations, pollen loads of visitors and a plant-level human-disturbance index.

> **Associated paper:** [Authors] ([Year]). Pollen limitation, flower visitors, and human disturbance in *Buddleja delavayi* (Scrophulariaceae). *[Journal]*. [DOI]
>
> **Archived version:** [Zenodo DOI, added after the first release]

## Quick start

1. Install R (tested with R 4.5.3) and the packages:

   ```r
   install.packages(c("readr", "dplyr", "tidyr", "ggplot2", "patchwork", "scales", "glmmTMB", "emmeans"))
   ```

2. Download or clone this repository and open `delavayi_pollination_pipeline.Rproj` in RStudio.
3. Run:

   ```r
   source("run_all.R")
   ```

The full run takes about two minutes and writes all results to `outputs/`.

## Repository structure

```
data/                          input data (four CSV files)
R/
  00_setup.R                   packages, settings, helper functions, figure style
  01_data.R                    reads and checks the data, builds the analysis tables
  02_pollen_limitation.R       fruit set, seed set, pollen limitation, years, repeated plants
  03_visitation.R              visitation rates, pollen limitation vs bee visitation
  04_pollen_loads.R            pollen loads, pollen-transport index, load composition
  05_disturbance.R             observer agreement, disturbance tests
  06_model_checks.R            residual checks of all models, sensitivity and robustness analyses
  07_figures.R                 figures
  08_tables.R                  Benjamini-Hochberg correction, identifier table, session info
run_all.R                      runs scripts 02-08 in order
outputs/                       results (created by the scripts)
```

Each script can also be run on its own once the earlier scripts have been run.

## Data

All files are comma-separated, UTF-8 encoded, with one header row. Missing values are empty cells.

### `pollination_experiment.csv` (one row per experimental unit)

| Column | Description |
|---|---|
| `year` | 2025 or 2026 |
| `plant_id` | permanent plant identifier (BD001–BD026), the same plant in both years |
| `field_label` | label used in the field that year (labels changed between years) |
| `pair` | number of the pair of units on that plant |
| `treatment` | `Control` (open pollination) or `Supplemental` (hand pollination added) |
| `flowers` | marked flowers in the unit |
| `fruits` | mature fruits |
| `ovules` | ovules in the dissected fruits (seeds + unfertilized ovules) |
| `seeds` | seeds in the dissected fruits |

### `visitation.csv` (one row per plant, year and visitor group)

| Column | Description |
|---|---|
| `year`, `plant_id`, `field_label` | as above |
| `inflorescences` | total inflorescences on the plant |
| `flowers_observed` | open flowers under observation |
| `visitor_group` | `Bees`, `Butterflies` or `Others` |
| `visits` | visits recorded |

### `disturbance.csv` (one row per plant and year)

| Column | Description |
|---|---|
| `year`, `plant_id`, `field_label` | as above |
| `road` | road proximity score, scored in the field |
| `herbivory_observer1`, `cutting_observer1` | field scores of the first observer |
| `herbivory_observer2`, `cutting_observer2` | scores from photographs by a second observer (used in the main analyses) |

All scores range from 0 to 2. Road proximity is 0 for more than 200 m, 1 for 50–200 m and 2 for less than 50 m from the nearest road. Herbivory is 0 for no evidence, 1 for some browsing marks or occasional dung and 2 for heavy browsing or active grazing. Cutting is 0 for no evidence, 1 for old weathered stumps and 2 for recent or repeated cutting.

### `pollen_loads.csv` (one row per insect)

| Column | Description |
|---|---|
| `insect_no` | unique number of the insect (1–83) |
| `year` | year of capture |
| `insect_label` | label used in the field |
| `plant_label` | label of the plant on which the insect was caught |
| `plant_id` | permanent identifier of that plant, empty if it was not a study plant |
| `insect_category` | category as recorded (Bee, Butterfly, Syrphid fly, Other) |
| `visitor_group` | `Bees`, `Butterflies` or `Others` |
| `b_delavayi` … `unidentified` | estimated number of pollen grains of each type (scaled from subsamples) |

Plant-years present in `disturbance.csv` but absent from `pollination_experiment.csv` are the three plants whose experiment was lost because they were cut. The data contain no plant coordinates.

## Outputs and the manuscript

| Manuscript item | Output |
|---|---|
| Table 1 | `outputs/tables/pollen_limitation_summary.csv` |
| Table 2 | `outputs/tables/pollen_load_summary.csv`, `pollen_load_estimates.csv`, `pollen_load_tests.csv`, `pollen_transport_index.csv`, `tests_BH.csv` |
| Fig. 1 | `outputs/figures/Fig1_visitor_rings_2025` (visitor rings, placed on the map drawn separately) |
| Figs. 2–5 | `outputs/figures/Fig2_pollen_limitation` … `Fig5_disturbance` |
| Figs. S1–S4 | `outputs/figures/FigS1_sampling_design` … `FigS4_pollen_load_composition` |
| Table S1 | `outputs/tables/plant_identifiers.csv` |
| Table S3 | `outputs/tables/interobserver_agreement.csv` |
| Table S4 | `outputs/model_checks/field_scores.csv`, `outputs/tables/disturbance_tests.csv` |
| Table S5 | `outputs/tables/pollen_composition.csv` |
| Table S6 | `outputs/model_checks/robustness.csv` |
| Appendix S1 | `outputs/model_checks/` (residual checks, residual plots, beta-binomial models, unit size) |

Figures are saved as PDF (editable text) and 600 dpi PNG, 180 mm wide. Fitted models are saved in `outputs/models/` together with the data they were fitted to. `outputs/session_info.txt` lists the R and package versions of the last run.

## Software

The results in the paper were produced with R 4.5.3, glmmTMB 1.1.15.2 and emmeans 2.0.3, together with readr, dplyr, tidyr, ggplot2, patchwork and scales. Bootstrap, permutation and simulation steps use a fixed random seed, so repeated runs give the same results.

## License

Code: [MIT License](LICENSE). Data: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Please cite the paper when using the data or code.

## Contact

[Name], [Institution], [e-mail]
