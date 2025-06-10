configfile: "config.yml"

SAMPLES = config["samples"]

rule all:
  input:
    "results/combined_samples_figures/filter_genes_dispersion_combined_samples.png",
    "results/combined_samples_figures/pca_variance_ratio_combined_samples.pdf",
    "results/combined_samples_figures/umap_neighbors_combined_samples.png",
    "results/combined_samples_figures/umap_leiden_combined_samples.pdf",
    "results/combined_samples_figures/umap_doublet_qc_combined_samples.pdf",
    "results/combined_samples_figures/umap_mt_counts_combined_samples.pdf",
    "results/combined_samples_figures/umap_resolutions_combined_samples.pdf",
    "results/combined_samples_figures/umap_resolutions_combined_samples.pdf",
    #"results/combined_samples_figures/dotplot_ranked_genes_combined_samples.png",
    "results/combined_samples_figures/umap_annotations_combined_samples.pdf",
    "results/combined_samples_figures/top10_degs_all_clusters_combined_samples.html",
    "results/combined_samples_figures/umap_celltypist_coarse_annotations.pdf",
    "results/combined_samples_figures/umap_celltypist_fine_annotations.pdf",


rule run_qc:
  input:
    h5ad = expand("results/{sample}_post_qc_filtered.h5ad", sample = SAMPLES)
  params:
    sample = "combined"
  output:
    highly_variable_genes = "results/combined_samples_figures/filter_genes_dispersion_combined_samples.png",
    pca_variance_ratio = "results/combined_samples_figures/pca_variance_ratio_combined_samples.pdf",
    neighbors = "results/combined_samples_figures/umap_neighbors_combined_samples.png",
    leiden = "results/combined_samples_figures/umap_leiden_combined_samples.pdf",
    reassess_qc_doublets = "results/combined_samples_figures/umap_doublet_qc_combined_samples.pdf",
    reassess_qc_mito = "results/combined_samples_figures/umap_mt_counts_combined_samples.pdf",
    umap_leiden_resolutions = "results/combined_samples_figures/umap_resolutions_combined_samples.pdf",
    umap_annotations = "results/combined_samples_figures/umap_annotations_combined_samples.pdf",
    #dotplot_ranked_genes = "results/combined_samples_figures/dotplot_ranked_genes_combined_samples.png",
    top10_degs = "results/combined_samples_figures/top10_degs_all_clusters_combined_samples.html",
    celltypist_coarse_umap = "results/combined_samples_figures/umap_celltypist_coarse_annotations.pdf",
    cell_typist_fine_umap = "results/combined_samples_figures/umap_celltypist_fine_annotations.pdf"
  shell:
    """
    python integrate_samples.py --input_h5ad_samples {input.h5ad}
    """

