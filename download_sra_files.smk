configfile: "config.yml"

SRR_IDS = config["srr_ids"]
SAMPLES = config["samples"]

rule all:
  input:
    expand("data/sra//{srr_id}_1.fastq.gz", srr_id = SRR_IDS),
    expand("data/sra//{srr_id}_2.fastq.gz", srr_id = SRR_IDS),
    expand("data/sra//{srr_id}_3.fastq.gz", srr_id = SRR_IDS),



rule fastq_dump:
  output:
    r1 = protected("data/sra//{srr_id}_1.fastq.gz"),
    r2 = protected("data/sra//{srr_id}_2.fastq.gz"),
    r3 = protected("data/sra//{srr_id}_3.fastq.gz")
  params:
    output_dir = "data/sra"
  shell:
    """
    fastq-dump {wildcards.srr_id} -x -O {params.output_dir}
    """

