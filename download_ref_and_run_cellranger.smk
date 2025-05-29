configfile: "config.yml"

SAMPLES = config["samples"]
SRR_IDS = config["srr_ids"]
FASTQ_DIR = config["fastq_dir"]
REFERENCE = config["reference"]

rule all:
    input:
        "data/cellranger/human_reference.tar.gz",
        "data/cellranger/cellranger-9.0.1.tar.gz",
        expand("data/cellranger/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", sample = SAMPLES)

rule download_reference:
    params:
        output_dir = "data/cellranger",
        url = "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
    output:
        output_file = "data/cellranger/human_reference.tar.gz"
    shell:
        """
        mkdir -p {params.output_dir}
        curl -o {output.output_file} {params.url}
        tar -xvzf {output.output_file} -C {params.output_dir}
        """

rule download_cell_ranger:
    params:
        output_dir = "data/cellranger",
        url = "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=1748433743&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=GgXUQvy4owG43kvJBcNwaOjvjGo3UcsKtfRf5tuXpnZVid5~w9Exks6pGP2TYC7SzFx76kTrfYgGx7TI-KA9AZ-XIUuuV4vvhBWe6PqwbSOwce7-fRpE0sYuQsSpfn7Ps-CxKnkIwxDvgCIcU9UfuJ8bQxEoGgAT6ewnyKpWJYDHZ~VhuqetkoRbpOCCKdzt37UiKE5Vhf9fX8CGw7Wq-MaF3nEbIV-V1-D24CtRhq6~W-hNmsnr~7BpSVTilYJgK5ywRKYKnJcJgnGBw7J18S8Dl4vk~Uqt6wQ1duKUkDFkxvE3Uty5I3imzxhqW0SODMp1Sa4YsKjaaETZnTCOBw__" 
    output:
        outfile = "data/cellranger/cellranger-9.0.1.tar.gz"
    shell:
        """
        mkdir -p {params.output_dir}
        curl -o {output.outfile} "{params.url}"
        tar -xvzf {output.outfile} -C {params.output_dir}
        """

rule run_cellranger_count:
  params:
    sample = "{sample}",
    fastqs = FASTQ_DIR,
    transcriptome = REFERENCE,
    cellranger_exe = "data/cellranger/cellranger-9.0.1/cellranger count"
  threads: 8
  output:
    output_file = "data/cellranger/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
  shell:
    """
    {params.cellranger_exe} --id="{params.sample}" --create-bam true \
    --transcriptome="{params.transcriptome}" \
    --fastqs="{params.fastqs}" \
    --sample="{params.sample}" \
    --localcores={threads}
    """
  
