configfile: "config.yml"

SAMPLES = config["samples"]

rule all:
  input:
    expand("results/{sample}_processed.h5ad", sample = SAMPLES),
    expand("results/figures/violin_{sample}.png", sample = SAMPLES),
    expand("results/figures/scatter_{sample}.png", sample = SAMPLES),
    expand("results/figures/filter_genes_dispersion_{sample}.png", sample = SAMPLES),
    expand("results/figures/pca_variance_ratio_{sample}.pdf", sample = SAMPLES),
    expand("results/figures/pca_{sample}.pdf", sample = SAMPLES),
    expand("results/figures/umap_neighbors_{sample}.png", sample = SAMPLES),
    expand("results/figures/umap_leiden_{sample}.pdf", sample = SAMPLES),
    expand("results/figures/umap_doublet_qc_{sample}.pdf", sample = SAMPLES),
    expand("results/figures/umap_mt_counts_{sample}.pdf", sample = SAMPLES),
    expand("results/{sample}_post_processed.h5ad", sample = SAMPLES),
    expand("results/figures/umap_resolutions_{sample}.pdf", sample = SAMPLES),
    expand("results/figures/dotplot_0.02_{sample}.pdf", sample = SAMPLES),
    expand("results/figures/dotplot_0.50_{sample}.pdf", sample = SAMPLES),
    expand("results/figures/umap_resolutions_{sample}.pdf", sample = SAMPLES),
    expand("results/figures/umap_annotations_{sample}.pdf", sample = SAMPLES),
    expand("results/figures/dotplot_0.02_{sample}.png", sample = SAMPLES),
    expand("results/figures/dotplot_0.50_{sample}.png", sample = SAMPLES),
    expand("results/figures/dotplot_ranked_genes_{sample}.png", sample = SAMPLES),
    expand("results/figures/top10_degs_all_clusters_{sample}.html", sample = SAMPLES),
rule run_scanpy:
  input:
    matrix_file = "{sample}_filtered_feature_bc_matrix.h5",
  output:
    h5ad = "results/{sample}_processed.h5ad"
  params:
    results_dir = "results",
    figures_dir = "results/figures"
  shell:
    """
    mkdir -p {params.results_dir}
    mkdir -p {params.figures_dir}
    python import_10x_data.py --feature_matrix_file {input.matrix_file} \
    --sample_id {wildcards.sample} --processed_h5 {output.h5ad}
    """

rule run_qc:
  input:
    h5ad = "results/{sample}_processed.h5ad"
  params:
    sample = "{sample}"
  output:
    violin = "results/figures/violin_{sample}.png",
    scatter_plot = "results/figures/scatter_{sample}.png",
    highly_variable_genes = "results/figures/filter_genes_dispersion_{sample}.png",
    pca_variance_ratio = "results/figures/pca_variance_ratio_{sample}.pdf",
    pca = "results/figures/pca_{sample}.pdf",
    neighbors = "results/figures/umap_neighbors_{sample}.png",
    leiden = "results/figures/umap_leiden_{sample}.pdf",
    reassess_qc_doublets = "results/figures/umap_doublet_qc_{sample}.pdf",
    reassess_qc_mito = "results/figures/umap_mt_counts_{sample}.pdf",
    umap_leiden_resolutions = "results/figures/umap_resolutions_{sample}.pdf",
    dotplot_05 = "results/figures/dotplot_0.50_{sample}.pdf",
    umap_annotations = "results/figures/umap_annotations_{sample}.pdf",
    dotplot_ranked_genes = "results/figures/dotplot_ranked_genes_{sample}.png",
    top10_degs = "results/figures/top10_degs_all_clusters_{sample}.html"
  shell:
    """
    python qc_scanpy.py --input_h5ad {input.h5ad} --sample_id {params.sample}
    """
