rule stage1_vep:
    conda:
        "../../envs/environment.yaml"
    input:
        vcf="data/processed/normalized.vcf.gz",
        vcf_index="data/processed/normalized.vcf.gz.csi",
        fasta=config["reference"]["fasta"],
        fasta_index=config["reference"]["fasta"] + ".fai"

    output:
        vcf="results/variants/stage1_minimal_vep.vcf.gz"

    log:
        "logs/stage1_vep.log"

    threads:
        config["stage1"]["vep"]["threads"]

    params:
        cache_dir=config["stage1"]["vep"]["cache_dir"],
        cache_version=config["stage1"]["vep"]["cache_version"],
        assembly=config["stage1"]["vep"]["assembly"],
        species=config["stage1"]["vep"]["species"],
        buffer_size=config["stage1"]["vep"]["buffer_size"]

    shell:
        r"""
        mkdir -p results/variants logs

        vep \
            --input_file {input.vcf} \
            --output_file {output.vcf} \
            --vcf \
            --compress_output bgzip \
            --cache \
            --offline \
            --dir_cache {params.cache_dir} \
            --cache_version {params.cache_version} \
            --species {params.species} \
            --assembly {params.assembly} \
            --fasta {input.fasta} \
            --fork {threads} \
            --buffer_size {params.buffer_size} \
            --symbol \
            --af_gnomadg \
            --pick \
            --fields "Allele,Consequence,SYMBOL,gnomADg_AF" \
            --no_stats \
            --force_overwrite \
            > {log} 2>&1
        """


rule stage1_route:
    conda:
        "../../envs/environment.yaml"
    input:
        vcf="results/variants/stage1_minimal_vep.vcf.gz",
        source="data/processed/normalized.vcf.gz",
        source_index="data/processed/normalized.vcf.gz.csi"

    output:
        vcf="results/variants/stage1_routed.vcf.gz"

    params:
        max_af=config["stage1"]["routing"]["max_gnomad_cache_af"]

    log:
        "logs/stage1_route.log"

    script:
        "../scripts/route_stage1.py"


rule index_stage1:
    conda:
        "../../envs/environment.yaml"
    input:
        vcf="results/variants/stage1_routed.vcf.gz"

    output:
        index="results/variants/stage1_routed.vcf.gz.csi"

    log:
        "logs/index_stage1.log"

    shell:
        r"""
        bcftools index \
            --force \
            --csi \
            {input.vcf} \
            > {log} 2>&1
        """
