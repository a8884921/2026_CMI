# 2026_CMI — figure-generation code

Source code for the figures of the CMI (2026) manuscript. Every figure is organised as a
directory of per-panel steps; each step is intended to be re-runnable on its own.

- Main figures: `Figure1` … `Figure6`
- Supplementary figures: `FigureS1` … `FigureS4`
- `privacy.R` — shared data-privacy guard sourced by every R panel.

## Repository layout

| Directory | Panels | Language | Step entry point |
|---|---|---|---|
| `Figure1/` | `02_PanelB` … `08_PanelH` | R | `run.sh` per panel |
| `Figure2/` | `01_PanelA1`, `02_PanelA2`, `03_PanelB` … `06_PanelE`, `08_Supplement` | R | `run.sh` per panel |
| `Figure3/` | `01_PanelA` … `08_PanelH` | R | `run.sh` per panel |
| `Figure4/` | `PanelA` … `PanelF` | R | `Rscript <panel>.R` (no wrapper) |
| `Figure5/` | `PanelA` … `PanelJ` | Python (matplotlib) | `run.sh` per panel |
| `Figure6/` | `PanelA` … `PanelG` | Python (matplotlib) | `run.sh` per panel |
| `FigureS1/` | `01_PanelA` … `04_PanelD` | R | `run.sh` per panel |
| `FigureS2/` | `PanelS4A` … `PanelS4D` | R | `Rscript <panel>.R` (no wrapper) |
| `FigureS3/` | `01_PanelA`, `02_PurityDoseResponse` | Python (matplotlib) | `run.sh` per panel |
| `FigureS4/` | `PanelA` … `PanelF` | Python (matplotlib) | `run.sh` per panel |

Layout inside a panel step:

- R panels: `<NN>_<panel>.R` + `run.sh` (plus `00_config.R` / extra scripts in `FigureS1`).
- Python panels: `render_panel.py` (entry point), `figure*.py` (panel drawing),
  `style*.py` (typography and colour theme), `data*.py` / `publication_data.py`
  (input reading and validation).

## Requirements

- **R** (≥ 4.1): `readr`, `dplyr`, `ggplot2`, `patchwork`, `tidyr`, `grid`, `pROC`,
  `lme4`, `circlize`, `ComplexHeatmap`.
  Set `RSCRIPT_BIN` if `Rscript` is not on `PATH`.
- **Python 3** (≥ 3.9): `numpy`, `pandas`, `matplotlib`, `scipy`, `pyarrow`
  (parquet inputs are read in the Python panels).
  Set `PYTHON_BIN` if `python3` is not on `PATH`.

## Running a panel

```bash
# R panel
bash Figure1/02_PanelB/run.sh

# Python panel
bash Figure6/PanelA/run.sh
```

Both wrappers `cd` into their own directory first, so they can be launched from anywhere.
Python panels honour `SCI_INPUT_DIR` (input tables) and `SCI_OUTPUT_DIR` (render output).

## Data and shared configuration

This repository contains the figure code only. Two pieces of the original analysis pipeline
are **not** bundled here and must be supplied before the panels run end to end:

1. `shared/figureN_config.R` — per-figure configuration helpers (exclusive to figures 1–4,
   S1 and S2). The R panels fail early with
   *“Required shared code is absent from the supplied archive”* when these are missing.
2. Step-local input tables — e.g. `Figure*/Panel*/inputs/plot_data/*.csv` read by the
   Python panels, and the panel CSV files read through each figure's config helper.

## Data privacy

`privacy.R` enforces the release data policy for every R panel:

- cohorts are referenced only through neutral aliases (`Cohort_*`, `Training`, `Test`, `All`);
- sample/patient identifiers must be pre-anonymised positive integers;
- free-text and direct-identifier columns (names, dates of birth, MRN, addresses, local file
  paths, notes … ) and non-English column names are rejected.

`publication_data.py` applies the analogous cohort-alias and identifier checks to the
Python panels.

## Citation

If you use this code, please cite the corresponding manuscript (CMI, 2026).
