configfile: "config.yml"
REFERENCE = config["reference"]


rule all:
  input:
    expand("data/human_reference/{reference}.zip", reference = REFERENCE),
    "data/cellranger/cellranger-9.0.1"


rule download_transcriptome:
  params:
    output_dir = "data/human_reference"
  output:
    output_file = "data/human_reference/{reference}.zip"
  shell:
    """
    mkdir -p {params.output_dir}
    datasets download genome accession {wildcards.reference} --filename {output.output_file} --include genome
    """

rule download_cell_ranger:
  params:
    output_dir = "data/cellranger",
    url = "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=1748238273&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=asgn8CF4jsasSZRWWfM3YQueRYAbP69VnfYTlQh-tk6wn1ZwzHJyV2C2R6iV4ERXjyoDmPAHIP9Klkwp2nd4hYwENrEIKiobdR~ATKcPsCRM3Dz05knra-7BbQDCJX8ZPQyYhqGQiqAmN8hKC9Z1YpiHF28EssRsQOsSy0-1mFwJ9AV3Hm3uhm32pyPKqe8Mzi36DdDvpf4wFFfdE08cR39zOYX~diA-NXOIBDUjvkjbJPg-JY3doNqQRYf8nLKg8F7h1y9Kdoh8AtaC2jY5sv8byuYPZ5S2jokA2rUtMeiySSS5MVL2nDMcv3rhWEdIqvNMfU~p2aYgsOgkcudyVQ__"
  output:
    outfile = "data/cellranger/cellranger-9.0.1.tar.gz"
  shell:
    """
    mkdir -p {params.output_dir}
    curl -o {output.outfile} "{params.url}"
    tar -xvzf {output.outfile} -C {params.output_dir}
    """
