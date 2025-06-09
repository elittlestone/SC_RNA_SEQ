config: "config.yml"

rule all:
  input:
    "results/combined_integrated.h5ad",
    "results/figures/umap_combined.png",
    "results/figures/leiden_combined.png"


rule integrate_samples:
  input:
    h5ads = expand("results/{sample}_processed.h5ad", sample = SAMPLES)
  output:
    integrated = "results/{integrated_sample_id}_integrated.h5ad",
    umap = "results/figures/umap_{integrated_sample_id}.png",
    leiden = "results/figures/leiden_{integrated_sample_id}.png"
  params:
    integrated_sample_id = config["integrated_sample_id"]
  shell:
    """
    python integrate_samples.py \
    --input_h5ads {input.h5ads} \
    --output_h5ad {output.integrated} \
    --sample_id {params.integrated_sample_id} \
    --umap_fig {output.umap} \
    --leiden_fig {output.leiden}
    """

