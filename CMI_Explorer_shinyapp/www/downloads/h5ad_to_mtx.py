#!/usr/bin/env python3
"""Convert an AnnData file (.h5ad) to the 10x Matrix Market ZIP that the Single-cell workspace reads.

    python3 h5ad_to_mtx.py data.h5ad S01.zip [--layer counts | --use-raw] [--where sample_id=S01 --where cell_type=B_cell]

Each --where KEY=VALUE[,VALUE...] keeps the cells whose obs column KEY has one of the values; several
--where filters are combined (one sample, B cells only). The older form --obs-key KEY --keep VALUE ...
still works.

Writes matrix.mtx.gz (genes x cells, integer counts), features.tsv.gz (gene ID, gene symbol,
"Gene Expression") and barcodes.tsv.gz into one ZIP. Requires: anndata, scipy, numpy.
Use raw UMI counts: adata.layers["counts"] or adata.raw.X when adata.X holds normalized values.
"""
import argparse
import gzip
import io
import sys
import zipfile

import anndata as ad
import numpy as np
import scipy.io
import scipy.sparse as sp


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("h5ad")
    p.add_argument("out_zip")
    p.add_argument("--layer", default=None, help="layer with raw counts (default: X, or raw.X if --use-raw)")
    p.add_argument("--use-raw", action="store_true", help="take counts from adata.raw")
    p.add_argument("--obs-key", default=None, help="obs column used to select cells (for example cell_type)")
    p.add_argument("--keep", nargs="*", default=None, help="values of --obs-key to keep (for example B_cell)")
    p.add_argument("--where", action="append", default=[], metavar="KEY=VALUE[,VALUE]",
                   help="keep cells whose obs column KEY has one of the values; repeat to combine filters")
    p.add_argument("--symbol-column", default=None, help="var column with gene symbols (default: var index)")
    a = p.parse_args()

    adata = ad.read_h5ad(a.h5ad)
    filters = [(a.obs_key, list(a.keep))] if a.obs_key and a.keep else []
    for w in a.where:
        if "=" not in w:
            sys.exit(f"--where needs KEY=VALUE, got: {w}")
        key, values = w.split("=", 1)
        filters.append((key, [v for v in values.split(",") if v]))
    for key, values in filters:
        if key not in adata.obs.columns:
            sys.exit(f"obs has no column '{key}'; available: {', '.join(map(str, adata.obs.columns))}")
        adata = adata[adata.obs[key].astype(str).isin(values)].copy()
    if adata.n_obs == 0:
        sys.exit("No cells remain after the filters.")
    if a.use_raw:
        X, var = adata.raw.X, adata.raw.var
    elif a.layer:
        X, var = adata.layers[a.layer], adata.var
    else:
        X, var = adata.X, adata.var
    X = sp.csr_matrix(X)
    if X.nnz and not np.allclose(X.data, np.round(X.data)):
        sys.exit("The selected matrix is not integer counts. Choose the raw-count layer with --layer or --use-raw.")
    X.data = np.round(X.data).astype(np.int64)
    symbols = var[a.symbol_column].astype(str).values if a.symbol_column else var.index.astype(str).values
    ids = var["gene_ids"].astype(str).values if "gene_ids" in var.columns else symbols

    def gz(text_or_bytes):
        b = text_or_bytes.encode() if isinstance(text_or_bytes, str) else text_or_bytes
        return gzip.compress(b, mtime=0)

    buf = io.BytesIO()
    scipy.io.mmwrite(buf, X.T.tocoo(), field="integer")  # genes x cells
    with zipfile.ZipFile(a.out_zip, "w", zipfile.ZIP_STORED) as z:
        z.writestr("matrix.mtx.gz", gz(buf.getvalue()))
        z.writestr("features.tsv.gz", gz("".join(f"{i}\t{s}\tGene Expression\n" for i, s in zip(ids, symbols))))
        z.writestr("barcodes.tsv.gz", gz("".join(f"{c}\n" for c in adata.obs_names.astype(str))))
    print(f"Wrote {a.out_zip}: {X.shape[1]} genes x {X.shape[0]} cells, {X.nnz} non-zero counts")


if __name__ == "__main__":
    main()
