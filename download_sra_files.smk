configfile: "config.yml"

SRR_IDS = config["srr_ids"]
SAMPLES = config["samples"]

rule all:
  input:
    expand("/Volumes/T7_Shield/temp/{srr_id}_1.fastq.gz", srr_id = SRR_IDS),
    expand("/Volumes/T7_Shield/temp/{srr_id}_2.fastq.gz", srr_id = SRR_IDS),
    expand("/Volumes/T7_Shield/temp/{srr_id}_3.fastq.gz", srr_id = SRR_IDS),



rule fastq_dump:
  output:
    r1 = protected("/Volumes/T7_Shield/temp/{srr_id}_1.fastq.gz"),
    r2 = protected("/Volumes/T7_Shield/temp/{srr_id}_2.fastq.gz"),
    r3 = protected("/Volumes/T7_Shield/temp/{srr_id}_3.fastq.gz")
  params:
    disk_limit = "200GB",
    temp_dir = "/Volumes/T7_Shield/temp",
    output_dir = "/Volumes/T7_Shield/temp"
  shell:
    """
    fastq-dump -x -O {params.output_dir} {wildcards.srr_id} --disk-limit {params.disk_limit}
    """

