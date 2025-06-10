import celltypist
from celltypist import models
import scanpy as sc
import scanpy.external as sce
import argparse
import scrublet as scr
import pandas as pd
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_h5ad_samples", nargs = "+", required=True)
    args = parser.parse_args()
    
    # Create dirs 
    outdir = "results/combined_samples_figures"
    os.makedirs(outdir, exist_ok=True)
    sc.settings.figdir = outdir

    adata_multiple_samples = args.input_h5ad_samples
    adata = combine_samples(adata_multiple_samples)
    # Post-processing steps go here
    normalization(adata)
    feature_selection(adata)
    dimensionality_reduction(adata)
    integration(adata)
    nearest_neighbor(adata)
    clustering(adata)
    reassess_qc(adata)
    annotation_and_identify_markers(adata)
    automated_annotation_celltypist(adata)

def combine_samples(h5ad_files):
    # Concatenate anndata objects from different samples into one
    adata_list = [sc.read_h5ad(file) for file in h5ad_files]
    sample_labels = [os.path.basename(file).split("_")[0] for file in h5ad_files]
    
    adata = sc.concat(
            adata_list,
            label = "sample",
            keys = sample_labels,
            index_unique = None
            )
    adata.obs_names_make_unique()
    return adata



def normalization(adata):
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    

def feature_selection(adata):
    sc.pp.highly_variable_genes(adata, n_top_genes = 2000, batch_key = "sample")
    sc.pl.highly_variable_genes(adata, save = f"_combined_samples.png")
    print(adata.var['highly_variable'].sum())
 

def dimensionality_reduction(adata):
    sc.tl.pca(adata)
    sc.pl.pca_variance_ratio(adata, n_pcs = 50, log = True, save = f"_combined_samples")
    
    # Could be useful if we had different samples in the same anndata object
    #sc.pl.pca(adata,
    #          color = ["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
    #          dimensions = [(0, 1), (2, 3), (0, 1), (2, 3)],
    #          ncols = 2,
    #          size = 2, 
    #          save = f"_combined_samples")
    

def integration(adata):
    sce.pp.harmony_integrate(adata, key = "sample")

def nearest_neighbor(adata):
    sc.pp.neighbors(adata, n_pcs = 10)
    sc.tl.umap(adata)
    sc.pl.umap(
            adata,
            color="sample",
            size=2,
            save=f"_neighbors_combined_samples.png"
            )
    

def clustering(adata):
    sc.tl.leiden(adata, flavor = "igraph", n_iterations = 2)
    sc.pl.umap(adata, color = ["leiden"], save = f"_leiden_combined_samples.pdf")


def reassess_qc(adata):
    sc.pl.umap(
            adata,
            color = ["leiden", "predicted_doublets", "doublet_scores"],
            wspace = 0.5, size = 3,
            save = f"_doublet_qc_combined_samples.pdf"
            )

    sc.pl.umap(adata,
               color = ["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
               wspace = 0.5, ncols = 2,
               save = f"_mt_counts_combined_samples.pdf")


def annotation_and_identify_markers(adata):

    # Decide which resolution to use for cell-type annotation 

    for res in [0.02, 0.5, 2.0]:
        sc.tl.leiden(adata,
                     key_added = f"leiden_res_{res:4.2f}", resolution = res, flavor = "igraph"
                     )
    sc.pl.umap(
            adata, 
            color = ["leiden_res_0.02", "leiden_res_0.50", "leiden_res_2.00"],
            legend_loc = "on data",
            save = f"_resolutions_combined_samples.pdf"
            )


    # Run DE analysis
    sc.tl.rank_genes_groups(adata,
                            groupby = "leiden_res_2.00",
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
    
    # Manual annotation dotplot
    #sc.pl.dotplot(adata, marker_genes_in_anndata, groupby = "leiden_res_0.50",
    #              standard_scale = "var",
    #              save = f"_combined_samples_0.50.pdf")


    #sc.tl.marker_gene_overlap(adata, marker_genes_in_anndata, key = "rank_genes_groups")


    #sc.pl.rank_genes_groups_dotplot(adata,
    #                                groupby = "leiden_res_0.50",
    #                                standard_scale = "var", n_genes = 5,
    #                                save = f"dotplot_ranked_genes_combined_samples.png")

    for category, genes in marker_genes_in_anndata.items():
        sc.tl.score_genes(adata, gene_list = genes, score_name = category)

    cluster_scores = (
            adata.obs.groupby("leiden_res_2.00")[list(marker_genes_in_anndata.keys())].mean()
            )

    top_cell_type_per_cluster = cluster_scores.idxmax(axis = 1)

    adata.obs["cell_type"] = adata.obs["leiden_res_2.00"].map(top_cell_type_per_cluster)

    sc.pl.umap(adata, color = "cell_type", legend_loc = "on data",
               save = f"_annotations_combined_samples.pdf")
    
    # Get top DE genes per cluster
    top_diff_exp_genes = pd.concat([
        sc.get.rank_genes_groups_df(adata, group = cl)
        .head(10)
        .assign(cluster = cl)
        for cl in adata.obs["leiden_res_2.00"].cat.categories
        ])

    top_diff_exp_genes.to_html(f"results/combined_samples_figures/top10_degs_all_clusters_combined_samples.html", index = False)

def automated_annotation_celltypist(adata):
    # Make copy of adata 
    adata_celltypist = adata.copy()
    # Convert .X to dense instead of sparse
    adata_celltypist.X = adata_celltypist.X.toarray()
    models.download_models(force_update = True,
                               model = ["Immune_All_Low.pkl", "Immune_All_High.pkl"]
                               )
    model_low = models.Model.load(model = "Immune_All_Low.pkl")
    model_high = models.Model.load(model = "Immune_All_High.pkl")
        
    # Get cell-type predictions
    predictions_high = celltypist.annotate(
            adata_celltypist, model = model_high, majority_voting = True
            )
    # Copy results back to original AnnData object
    predictions_high_adata = predictions_high.to_adata()
        
    # Add annotations for high and low level back to original adata object
    adata.obs["celltypist_cell_label_coarse"] = predictions_high_adata.obs.loc[
            adata.obs.index, "majority_voting"
            ]
    adata.obs["celltypist_conf_score_coarse"] = predictions_high_adata.obs.loc[
            adata.obs.index, "conf_score"
            ]
                        
    predictions_low = celltypist.annotate(
            adata_celltypist, model = model_low, majority_voting = True
            )
    predictions_low_adata = predictions_low.to_adata()
    
        
    adata.obs["celltypist_cell_label_fine"] = predictions_low_adata.obs.loc[
            adata.obs.index, "majority_voting"
            ]
    adata.obs["celltypist_conf_score_fine"] = predictions_low_adata.obs.loc[
            adata.obs.index, "conf_score"
            ]
        
    # Now plot !
    sc.pl.umap(
            adata,
            color = ["celltypist_cell_label_coarse", "celltypist_conf_score_coarse"],
            frameon = False,
            sort_order = False,
            wspace = 1,
            legend_loc = "on data",
            save = f"_celltypist_coarse_annotations.pdf"
            )
        
    sc.pl.umap(
            adata,
            color = ["celltypist_cell_label_fine", "celltypist_conf_score_fine"],
            frameon = False,
            sort_order = False,
            wspace = 1,
            save = f"_celltypist_fine_annotations.pdf"
            )


if __name__ == "__main__":
    main()

