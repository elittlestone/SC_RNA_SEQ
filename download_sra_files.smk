configfile: "config.yml"

SRR_IDS = config["srr_ids"]

rule all:
  input:
    expand("data/sra/{srr_id}", srr_id = SRR_IDS)


rule download_reads:
  params:
    output_dir = "data/sra",
    max_size = "25G"
  output:
    sra_files = protected(directory("data/sra/{srr_id}"))
  shell:
    """fasterq-dump -p -v --split-3 -O {params.output_dir} {wildcards.srr_id} 
    --disk-limit {params.max_size} -t /tmp"""

