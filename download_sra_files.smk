configfile: "config.yml"

SRR_IDS = config["srr_ids"]
SAMPLES = config["samples"]

PAIRS = list(zip(SAMPLES, SRR_IDS))

rule all:
  input:
    expand("data/sra/{sample}_S1_R1_001.fastq", sample = SAMPLES)
    expand("data/sra/{sample}_S1_R2_001.fastq", sample = SAMPLES)

rule fasterq_dump:
  output:
    r1 = protected("data/sra/{srr_id}_1.fastq"),
    r2 = protected("data/sra/{srr_id}_2.fastq"),
    r3 = protected("data/sra/{srr_id}_3.fastq")
  params:
    output_dir = "data/sra"
  shell:
    """
    fasterq-dump {wildcards.srr_id} -x -O {params.output_dir}
    """

