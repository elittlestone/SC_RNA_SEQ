import scanpy as sc
import argparse
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_h5ad", required=True)
    parser.add_argument("--sample_id", required=True)
    args = parser.parse_args()
    
    # Create dirs 
    outdir = "results/figures"
    os.makedirs(outdir, exist_ok=True)
    sc.settings.figdir = outdir

    adata = sc.read_h5ad(args.input_h5ad)
    sample_id = args.sample_id
    # Post-processing steps go here
    adata = filter_cells(adata, sample_id)
    adata = doublet_detection(adata)
    adata = normalization(adata)
    adata = feature_selection(adata, sample_id)
    adata = dimensionality_reduction(adata, sample_id)



def filter_cells(adata, sample_id):

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
        save=f"_{sample_id}.png"
    )

    sc.pl.scatter(
        adata,
        x="total_counts",
        y="n_genes_by_counts",
        color="pct_counts_mt",
        save=f"_{sample_id}.png"
    )
    
    sc.pp.filter_cells(adata, min_genes = 100)
    sc.pp.filter_genes(adata, min_cells = 3)

    return adata

def doublet_detection(adata):
    sc.pp.scrublet(adata, batch_key="sample")
    return adata 

def normalization(adata):
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    return adata

def feature_selection(adata, sample_id):
    sc.pp.highly_variable_genes(adata, n_top_genes = 2000, batch_key = "sample")
    sc.pl.highly_variable_genes(adata, save = f"_{sample_id}.png")
    return adata

def dimensionality_reduction(adata, sample_id):
    sc.tl.pca(adata)
    sc.pl.pca_variance_ratio(adata, n_pcs = 50, log = True, save = f"_{sample_id}")
    sc.pl.pca(adata,
              color = ["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
              dimensions = [(0, 1), (2, 3), (0, 1), (2, 3)],
              ncols = 2,
              size = 2, 
              save = f"_{sample_id}")
    return adata
if __name__ == "__main__":
    main()

