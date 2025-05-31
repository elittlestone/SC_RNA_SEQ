configfile: "config.yml"

SAMPLES = config["samples"]

rule all:
  input:
    expand("results/{sample}_processed.h5ad", sample = SAMPLES),
    expand("results/figures/violin_{sample}.png", sample = SAMPLES),
    expand("results/figures/scatter_{sample}.png", sample = SAMPLES),
    expand("results/figures/filter_genes_dispersion_{sample}.png", sample = SAMPLES),
    expand("results/figures/pca_variance_ratio_{sample}.pdf", sample = SAMPLES),
    expand("results/figures/pca_{sample}.pdf", sample = SAMPLES)

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
    pca = "results/figures/pca_{sample}.pdf"
  shell:
    """
    python qc_scanpy.py --input_h5ad {input.h5ad} --sample_id {params.sample}
    """
