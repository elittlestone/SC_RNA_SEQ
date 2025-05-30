configfile: "config.yml"

SAMPLES = config["samples"]

rule all:
  input:
    expand("results/{sample}_processed.h5ad", sample = SAMPLES),
    expand("results/qc_figures/violin_{sample}.png", sample = SAMPLES),
    expand("results/qc_figures/scatter_{sample}.png", sample = SAMPLES)
        

rule run_scanpy:
  input:
    matrix_file = "{sample}_filtered_feature_bc_matrix.h5",
  output:
    h5ad = "results/{sample}_processed.h5ad"
  params:
    results_dir = "results",
    qc_figures_dir = "results/qc_figures"
  shell:
    """
    mkdir -p {params.results_dir}

    python run_scanpy.py --feature_matrix_file {input.matrix_file} \
    --sample_id {wildcards.sample} --processed_h5 {output.h5ad}
    """

rule run_qc:
  input:
    h5ad = "results/{sample}_processed.h5ad"
  params:
    sample = "{sample}",
    qc_dir = "results/qc_figures"
  output:
    violin = "results/qc_figures/violin_{sample}.png",
    scatter_plot = "results/qc_figures/scatter_{sample}.png"
  shell:
    """
    python qc_scanpy.py --input_h5ad {input.h5ad} --sample_id {params.sample}
    """
