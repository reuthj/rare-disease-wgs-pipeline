# =========================================================
# Reference validation
# =========================================================


rule index_reference:
    conda:
        "../../envs/environment.yaml"
    input:
        fasta=config["reference"]["fasta"]

    output:
        fai=config["reference"]["fasta"] + ".fai"

    log:
        "logs/index_reference.log"

    shell:
        r"""
        set -euo pipefail

        samtools faidx {input.fasta} \
            2> {log}
        """


rule validate_reference:
    conda:
        "../../envs/environment.yaml"
    input:
        header="results/qc/input_header.txt",
        fai=config["reference"]["fasta"] + ".fai"

    output:
        report="results/qc/reference_validation.tsv"

    script:
        "../scripts/validate_reference.py"
