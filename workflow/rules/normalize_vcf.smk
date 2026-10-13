# =========================================================
# VCF normalization
# =========================================================


rule normalize_vcf:
    conda:
        "../../envs/environment.yaml"
    input:
        vcf=config["input_vcf"],
        ref=config["reference"]["fasta"],
        fai=config["reference"]["fasta"] + ".fai"

    output:
        vcf="data/processed/normalized.vcf.gz"

    log:
        "logs/normalize_vcf.log"

    threads:
        2

    shell:
        r"""
        set -euo pipefail

        bcftools norm \
            -f {input.ref} \
            -c e \
            -m -any \
            --old-rec-tag ORIG \
            --threads {threads} \
            -Oz \
            -o {output.vcf} \
            {input.vcf} \
            2> {log}
        """


rule index_normalized_vcf:
    conda:
        "../../envs/environment.yaml"
    input:
        vcf="data/processed/normalized.vcf.gz"

    output:
        index="data/processed/normalized.vcf.gz.csi"

    log:
        "logs/index_normalized_vcf.log"

    shell:
        r"""
        set -euo pipefail

        bcftools index \
            --csi \
            --force \
            {input.vcf} \
            2> {log}
        """


rule normalized_stats:
    conda:
        "../../envs/environment.yaml"
    input:
        vcf="data/processed/normalized.vcf.gz",
        index="data/processed/normalized.vcf.gz.csi"

    output:
        stats="results/qc/normalized_stats.txt"

    log:
        "logs/normalized_stats.log"

    shell:
        r"""
        set -euo pipefail

        bcftools stats {input.vcf} \
            > {output.stats} \
            2> {log}
        """
