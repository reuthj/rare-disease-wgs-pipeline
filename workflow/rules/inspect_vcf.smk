rule index_vcf:
    conda:
        "../../envs/environment.yaml"
    input:
        vcf=config["input_vcf"]

    output:
        index=config["input_vcf"] + ".csi"

    log:
        "logs/index_vcf.log"

    threads:
        2

    shell:
        r"""
        set -euo pipefail

        gzip -t {input.vcf} 2> {log}

        bcftools index \
            --csi \
            --force \
            --threads {threads} \
            {input.vcf} \
            2>> {log}
        """


rule inspect_vcf:
    conda:
        "../../envs/environment.yaml"
    input:
        vcf=config["input_vcf"],
        index=config["input_vcf"] + ".csi"

    output:
        header="results/qc/input_header.txt",
        samples="results/qc/input_samples.txt",
        stats="results/qc/input_bcftools_stats.txt"

    log:
        "logs/input_inspect_vcf.log"

    shell:
        r"""
        set -euo pipefail

        bcftools view -h {input.vcf} \
            > {output.header} \
            2> {log}

        bcftools query -l {input.vcf} \
            > {output.samples} \
            2>> {log}

        bcftools stats {input.vcf} \
            > {output.stats} \
            2>> {log}
        """


rule summarize_vcf:
    conda:
        "../../envs/environment.yaml"
    input:
        vcf=config["input_vcf"],
        index=config["input_vcf"] + ".csi",
        header="results/qc/input_header.txt",
        samples="results/qc/input_samples.txt",
        stats="results/qc/input_bcftools_stats.txt"

    output:
        summary="results/qc/input_summary.tsv",
        fields="results/qc/input_header_fields.tsv"

    script:
        "../scripts/summarize_vcf.py"
