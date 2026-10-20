# Example data

Example inputs for each analysis and the result bundles this version computes for them. The files
are simulated or are software test fixtures; they show how the software works and are not biology.

| Analysis | File | Contents |
|---|---|---|
| Single-cell | `single_cell/simulated_b_cells.zip` | 10x Matrix Market ZIP, 36,601 genes × 900 simulated cells. Shown by "Try the example". |
| | `single_cell/counts_matrix_market.zip` | 10x Matrix Market ZIP, 36,601 genes × 96 synthetic cells (software test). |
| | `single_cell/counts.csv` | The same test cells as a genes × cells CSV. |
| Spatial | `spatial/spot_weights_demo.csv` | Spot table: spot_id, sample_id, x, y and 12 class weights; 2,171 simulated spots. Shown by "Try the example". |
| | `../www/downloads/rctd_test_data.zip` | Small reference (300 cells) with labels and a 36-spot section, for the RCTD script. |
| | `spatial/rctd_bundle_demo.zip` | Five-file RCTD bundle (local edition). |
| Bulk | `bulk/expression.csv` | Genes × samples CSV, 5,000 genes × 24 samples. Shown by "Try the example". |

`precomputed/<analysis>/` holds the result bundle of each example (made by
`tools/precompute_examples.R`), which "Try the example" displays at once.

## Expected results

- **Single-cell** (simulated B cells): 798 of 900 cells (88.7%) at or above the CMI threshold; 263 Mal-1,
  343 Mal-2 and 192 Mal-3 cells.
- **Spatial** (simulated section): 343 Mal-1, 230 Mal-2 and 161 Mal-3 spots; 1,437 non-malignant spots.
- **Bulk**: all 100 gene pairs computable; 12 CMI-high and 12 CMI-low samples (median CMI-riskscore 36.40).
