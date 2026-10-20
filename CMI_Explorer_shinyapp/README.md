# CMI Explorer 2.0.0

Computes the clonal malignancy index (CMI) and its grades for B-cell lymphoma data, from fixed model
parameters. The site has three pages: **Analysis**, **Help** and **About**.

| Analysis | Input | Output |
|---|---|---|
| Single-cell | Raw UMI counts of B cells from one sample (10x Matrix Market ZIP or CSV/TSV) | CMI and CMI grade (Below threshold, Mal-1, Mal-2, Mal-3) of every cell |
| Spatial | Spot table with deconvolution weights for Mal-1, Mal-2, Mal-3 and other classes (one section); `www/downloads/rctd_mal_states.R` makes it with RCTD | CMI grade of every spot (Mal-1, Mal-2, Mal-3 or non-malignant) |
| Bulk | Gene × sample expression matrix | CMI-riskscore of every sample and CMI-high/low grade at the cohort median |

| Grade | CMI |
|---|---|
| Below threshold | < 0.9755 |
| Mal-1 | 0.9755–0.98926 |
| Mal-2 | 0.98926–0.99620 |
| Mal-3 | > 0.99620 |

Every run returns a result bundle: tables, figures (SVG, PDF, 600 dpi PNG), quality-control files, a
README and `provenance.json`. Research use only; not for diagnosis or treatment decisions.

## Folder

```
app.R         interface (shiny + bslib); analyses run in a separate R process
R/            workspaces (Analysis page), Help and About pages, upload checks, figures, job control
engine/       analysis engines and the worker script (engine/worker.R)
ref/          model parameters and the reference files the calculations need
www/          styles, script, fonts, downloads/ (RCTD script and test data, cut points, h5ad_to_mtx.py)
examples/     example inputs and their precomputed results
tools/        deployment, package installation, example builder, integrity check
tests/        verification against the original implementation
run_local.R   local edition (files up to 4 GB, no cell or spot limits, RCTD available)
```

## Run

R 4.3 or later.

```sh
Rscript tools/install_packages.R      # once; add --rctd for RCTD in the local edition
Rscript -e 'shiny::runApp(".")'       # hosted edition, with the shinyapps.io limits
Rscript run_local.R                   # local edition
```

Many samples without the browser (one call per sample; same outputs as the site):

```sh
CMI_PROFILE=local Rscript engine/worker.R --root . --module single_cell \
  --input sample1.zip --output results/sample1 --state results/sample1/state.json
```

## Deploy to shinyapps.io

```r
install.packages("rsconnect")
rsconnect::setAccountInfo(name = "<account>", token = "<token>", secret = "<secret>")
```
```sh
Rscript tools/deploy_shinyapps.R --dry-run             # every check and the dependency scan; nothing is uploaded
Rscript tools/deploy_shinyapps.R <account> cmi-explorer
```

The script checks the model files against `ref/asset_manifest.json`, the precomputed examples and
`tools/check_integrity.R`, uploads only what the site needs (`app.R`, `R/`, `engine/`, `ref/`, `www/`,
`examples/`) and sets the instance size to Large (1 GB). Set **Max Worker Processes = 1** in the
dashboard. Details: `docs/DEPLOYMENT_zh.md`.

## Maintenance

| After changing | Run |
|---|---|
| `ref/` or `engine/constants.R` | `Rscript tools/export_gene_sets.R .` |
| `engine/`, `R/figures.R` or `examples/` | `Rscript tools/precompute_examples.R .` |
| anything numerical | `Rscript tests/test_equivalence.R .` |
| `www/downloads/rctd_mal_states.R` | `Rscript tests/test_rctd_script.R .` |
| anything | `Rscript tools/check_integrity.R .` |

After an approved change to a locked file, run the tests and then
`Rscript tools/check_integrity.R . --write-locks`.

Fonts: IBM Plex and Liberation Sans, both under the SIL Open Font License (`www/fonts/`, `engine/fonts/`).
