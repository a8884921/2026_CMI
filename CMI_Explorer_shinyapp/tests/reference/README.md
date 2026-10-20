# Reference outputs of the original release (v1.2.0)

The three files in this folder are the `results.csv` tables of the v1.2.0 result bundles in the
original package (`reports/local_e2e/*_result_bundle.zip`), copied byte for byte. They were produced
by the v1.2.0 local edition on 2026-10-16 from the example inputs that
are shipped here unchanged (SHA-256 checked by `tests/test_equivalence.R`).

| File | v1.2.0 bundle (SHA-256) | Input (SHA-256 of the input file) | Run |
|---|---|---|---|
| `v1.2.0_single_cell_results.csv` | `single_cell_result_bundle.zip` (c43182f5c4271f38...) | `examples/single_cell/counts.csv` (4beb978da1fcae9d...) | 2026-10-16 13:38-13:39 UTC |
| `v1.2.0_bulk_results.csv` | `bulk_result_bundle.zip` (80088d57586c78f2...) | `examples/bulk/expression.csv` (00a3905f1c2c48b7...) | 2026-10-16 13:39 UTC |
| `v1.2.0_spatial_rctd_results.csv` | `spatial_result_bundle.zip` (979a16fa72184854...) | `examples/spatial/rctd_bundle_demo.zip` (f82ef5c6f589cee9...) | 2026-10-16 13:39 UTC |

Environment recorded in the bundles' `provenance.json` and `R_session_info.txt`:

- R 4.4.1 (aarch64-apple-darwin20, macOS 15.1); app_version "1.2.0-local"
- Single cell: Seurat 5.4.0, SeuratObject 5.3.0, UCell 2.10.1, caret 7.0-1, gbm 2.2.3, qlcMatrix 0.9.9,
  org.Hs.eg.db 3.20.0, AnnotationDbi 1.68.0, h2o 3.46.0.2, Matrix 1.7-4; CNV with Python 3.10.18,
  numpy 1.26.4, scipy 1.15.2, anndata 0.11.4, scanpy 1.11.4, infercnvpy 0.6.1
- Bulk: randomForestSRC 3.5.0, HGNChelper 0.8.15, data.table 1.18.2.1, matrixStats 1.5.0
- Spatial: spacexr 2.2.1 (RCTD, full mode), Matrix 1.7-4

The same folder of the original package also contains `single_cell_reference_results.csv`, whose
`cnv_score` column differs from the bundle by up to 1e-4 (CMI by up to 0.0012, the same malignant
calls). It is not a result bundle of the release and is not used here; the bundle above is the output
of the documented infercnvpy path.

`tests/test_equivalence.R` compares the engines of this version with these files (tolerances: 1e-12
for the five dimensions, 2e-9 for the CMI because H2O stores its trees in single precision, 1e-10 for
the CMI-riskscore and the RCTD weights) and writes `tests/results/equivalence_report.json`.
