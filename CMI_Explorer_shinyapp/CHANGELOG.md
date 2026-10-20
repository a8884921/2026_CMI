# Changelog

## 2.0.0 (2026-10-18)

- Three pages: Analysis (Single-cell, Spatial and Bulk, switched with buttons), Help and About.
- Single-cell: the CMI of each cell and its grade (Below threshold, Mal-1, Mal-2, Mal-3).
- Spatial: the CMI grade of each spot (Mal-1, Mal-2, Mal-3 or non-malignant) from deconvolution weights.
- Bulk: the CMI-riskscore of each sample and the CMI-high/low grade at the cohort median.
- The package holds only the model parameters the calculations need; all calculations run in R.
- Uploads are checked on arrival; each analysis runs in its own R process with a first-come,
  first-served queue; results come as tables, figures (SVG, PDF, 600 dpi PNG) and a result bundle.
- Local edition (`run_local.R`): files up to 4 GB, no cell, spot or value limits, RCTD deconvolution.

## 1.2.0

First release (R, Python and Java engines).
