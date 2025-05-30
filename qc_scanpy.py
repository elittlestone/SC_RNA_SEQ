import scanpy as sc
import argparse
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_h5ad", required=True)
    parser.add_argument("--sample_id", required=True)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input_h5ad)

    # Set output dir
    outdir = "results/qc_figures"
    os.makedirs(outdir, exist_ok=True)
    sc.settings.figdir = outdir

    # Calculate QC metrics
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

    # Save plots
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        save=f"_{args.sample_id}.png"
    )

    sc.pl.scatter(
        adata,
        x="total_counts",
        y="n_genes_by_counts",
        color="pct_counts_mt",
        save=f"_{args.sample_id}.png"
    )

if __name__ == "__main__":
    main()

