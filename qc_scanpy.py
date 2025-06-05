import scrublet as scr
import scanpy as sc
import argparse
import pandas as pd
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
    filter_cells(adata, sample_id)
    doublet_detection(adata)
    normalization(adata)
    feature_selection(adata, sample_id)
    dimensionality_reduction(adata, sample_id)
    nearest_neighbor(adata, sample_id)
    clustering(adata, sample_id)
    reassess_qc(adata, sample_id)
    annotation_and_identify_markers(adata, sample_id)

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

    

def doublet_detection(adata):
    scrub = scr.Scrublet(adata.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    adata.obs["doublet_scores"] = doublet_scores
    adata.obs["predicted_doublets"] = predicted_doublets

def normalization(adata):
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    

def feature_selection(adata, sample_id):
    sc.pp.highly_variable_genes(adata, n_top_genes = 2000, batch_key = "sample")
    sc.pl.highly_variable_genes(adata, save = f"_{sample_id}.png")
    

def dimensionality_reduction(adata, sample_id):
    sc.tl.pca(adata)
    sc.pl.pca_variance_ratio(adata, n_pcs = 50, log = True, save = f"_{sample_id}")
    #sc.pl.pca(adata,
    #          color = ["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
    #          dimensions = [(0, 1), (2, 3), (0, 1), (2, 3)],
    #          ncols = 2,
    #          size = 2, 
    #          save = f"_{sample_id}")
    

    

def nearest_neighbor(adata, sample_id):
    sc.pp.neighbors(adata, n_pcs = 10)
    sc.tl.umap(adata)
    sc.pl.umap(
            adata,
            color="sample",
            size=2,
            save=f"_neighbors_{sample_id}.png"
            )
    

def clustering(adata, sample_id):
    sc.tl.leiden(adata, flavor = "igraph", n_iterations = 2)
    sc.pl.umap(adata, color = ["leiden"], save = f"_leiden_{sample_id}")


def reassess_qc(adata, sample_id):
    sc.pl.umap(
            adata,
            color = ["leiden", "predicted_doublets", "doublet_scores"],
            wspace = 0.5, size = 3,
            save = f"_doublet_qc_{sample_id}.pdf"
            )

    sc.pl.umap(adata,
               color = ["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
               wspace = 0.5, ncols = 2,
               save = f"_mt_counts_{sample_id}.pdf")


def annotation_and_identify_markers(adata, sample_id):
    for res in [0.02, 0.5, 2.0]:
        sc.tl.leiden(adata,
                     key_added = f"leiden_res_{res:4.2f}", resolution = res, flavor = "igraph"
                     )
    sc.pl.umap(
            adata, 
            color = ["leiden_res_0.02", "leiden_res_0.50", "leiden_res_2.00"],
            legend_loc = "on data",
            save = f"_resolutions_{sample_id}.pdf"
            )


    # Run DE analysis
    sc.tl.rank_genes_groups(adata,
                            groupby = "leiden_res_0.50",
                            method = "wilcoxon")
    
    # Dotplot and manual annotation

    marker_genes_all = {
    "CD14+ Mono" : ["FCN1", "CD14"], 
    "CD16+ Mono": ["TCF7L2", "FCGR3A", "LYN"],
    # Note: DMXL2 should be negative
    "cDC2": ["CST3", "COTL1", "LYZ", "DMXL2", "CLEC10A", "FCER1A"],
    "Erythroblast": ["MKI67", "HBA1", "HBB"],
    # Note HBM and GYPA are negative markers
    "Proerythroblast": ["CDK6", "SYNGR1", "HBM", "GYPA"],
    "NK": ["GNLY", "NKG7", "CD247", "FCER1G", "TYROBP", "KLRG1", "FCGR3A"],
    "ILC": ["ID2", "PLCG2", "GNLY", "SYNE1"],
    "Naive CD20+ B": ["MS4A1", "IL4R", "IGHD", "FCRL1", "IGHM"],
    # Note IGHD and IGHM are negative markers
    "B cells" : ["MS4A1", "ITGB1", "COL4A4", "PRDM1", "IRF4", "PAX5", "BCL11A", "BLK", "IGHD", "IGHM"],
    "Plasma cells": ["MZB1", "HSP90B1", "FNDC3B", "PRDM1", "IGKC", "JCHAIN"],
    # Note PAX5 is a negative marker
    "Plasmablast": ["XBP1", "PRDM1", "PAX5"],
    "CD4+ T": ["CD4", "IL7R", "TRBC2"],
    "CD8+ T": ["CD8A", "CD8B", "GZMK", "GZMA", "CCL5", "GZMB", "GZMH", "GZMA"],
    "T naive": ["LEF1", "CCR7", "TCF7"],
    "pDC": ["GZMB", "IL3RA", "COBLL1", "TCF4"]}

    marker_genes_in_anndata = {cell_type: [gene for gene in genes if gene in adata.var_names]
                               for cell_type, genes in marker_genes_all.items()}
    

    sc.pl.dotplot(adata, marker_genes_in_anndata, groupby = "leiden_res_0.50",
                  standard_scale = "var",
                  save = f"_0.50_{sample_id}.pdf")


    sc.tl.marker_gene_overlap(adata, marker_genes_in_anndata, key = "rank_genes_groups")


    for category, genes in marker_genes_in_anndata.items():
        sc.tl.score_genes(adata, gene_list = genes, score_name = category)


    sc.pl.umap(adata, color = list(marker_genes_in_anndata.keys()),
               save = f"_annotations_{sample_id}.pdf")


    sc.pl.rank_genes_groups_dotplot(adata,
                                    groupby = "leiden_res_0.50",
                                    standard_scale = "var", n_genes = 5,
                                    save = f"ranked_genes_{sample_id}.png")

    
    # Get top DE genes per cluster
    top_diff_exp_genes = pd.concat([
        sc.get.rank_genes_groups_df(adata, group = cl)
        .head(10)
        .assign(cluster = cl)
        for cl in adata.obs["leiden_res_0.50"].cat.categories
        ])

    top_diff_exp_genes.to_html(f"results/figures/top10_degs_all_clusters_{sample_id}.html", index = False)

if __name__ == "__main__":
    main()

