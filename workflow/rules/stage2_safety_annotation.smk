rule stage2_vep_context:
    conda:
        "../../envs/environment.yaml"
    input:
        vcf="results/variants/stage1_routed.vcf.gz",
        vcf_index="results/variants/stage1_routed.vcf.gz.csi",
        fasta=config["reference"]["fasta"],
        fasta_index=config["reference"]["fasta"] + ".fai",
        ccre=config["stage2"]["encode"]["ccre_bed"],
        ccre_index=config["stage2"]["encode"]["ccre_bed"] + ".tbi"

    output:
        temp("results/intermediate/stage2/context.vcf.gz")

    log:
        "logs/stage2_vep.log"

    threads:
        config["stage1"]["vep"]["threads"]

    params:
        cache_dir=config["stage1"]["vep"]["cache_dir"],
        cache_version=config["stage1"]["vep"]["cache_version"],
        assembly=config["stage1"]["vep"]["assembly"],
        species=config["stage1"]["vep"]["species"],
        buffer_size=config["stage1"]["vep"]["buffer_size"],
        pick_order=config["stage2"]["vep"]["pick_order"]

    shell:
        r"""
        mkdir -p results/intermediate/stage2 logs

        vep \
            --input_file {input.vcf} \
            --output_file {output} \
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
            --biotype \
            --canonical \
            --mane \
            --numbers \
            --transcript_version \
            --regulatory \
            --flag_pick_allele_gene \
            --pick_order {params.pick_order} \
            --custom file={input.ccre},short_name=ENCODE_cCRE,format=bed,type=overlap,coords=0 \
            --fields "Allele,Consequence,SYMBOL,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,MANE_SELECT,PICK,ENCODE_cCRE" \
            --no_stats \
            --force_overwrite \
            > {log} 2>&1
        """


rule stage2_add_clinvar:
    conda:
        "../../envs/environment.yaml"
    input:
        vcf="results/intermediate/stage2/context.vcf.gz",
        clinvar=config["stage2"]["clinvar"]["vcf"],
        clinvar_index=config["stage2"]["clinvar"]["vcf"] + ".tbi"

    output:
        temp("results/intermediate/stage2/context_clinvar.vcf.gz")

    log:
        "logs/stage2_clinvar.log"

    shell:
        r"""
        mkdir -p results/intermediate/stage2

        TMP=$(mktemp -d)

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/CLNSIG\n' \
            {input.clinvar} \
        | bgzip -c > "$TMP/clinvar.tsv.gz"

        tabix -f -s 1 -b 2 -e 2 "$TMP/clinvar.tsv.gz"

        printf '%s\n' \
          '##INFO=<ID=CLINVAR_ID,Number=1,Type=String,Description="ClinVar Variation ID">' \
          '##INFO=<ID=CLINVAR_CLNSIG,Number=1,Type=String,Description="ClinVar clinical significance">' \
          > "$TMP/header.txt"

        bcftools annotate \
            -a "$TMP/clinvar.tsv.gz" \
            -h "$TMP/header.txt" \
            -c CHROM,POS,REF,ALT,INFO/CLINVAR_ID,INFO/CLINVAR_CLNSIG \
            -Oz \
            -o {output} \
            {input.vcf} \
            > {log} 2>&1

        rm -rf "$TMP"
        """

rule stage2_build_safety:
    conda:
        "../../envs/environment.yaml"
    input:
        vcf="results/intermediate/stage2/context_clinvar.vcf.gz",
        source="results/variants/stage1_routed.vcf.gz",
        source_index="results/variants/stage1_routed.vcf.gz.csi"

    output:
        parquet="results/tables/stage2_safety.parquet"

    script:
        "../scripts/build_stage2_safety.py"

rule complete_stage2_safety:
    conda:
        "../../envs/environment.yaml"
    input:
        parquet="results/tables/stage2_safety.parquet",
        vcf="data/processed/normalized.vcf.gz"
    output:
        "results/tables/stage2_safety_complete.parquet"
    script:
        "../scripts/add_stage2_genotype.py"

rule stage2_variant_keys:
    conda:
        "../../envs/environment.yaml"
    input:
        vcf="results/variants/stage1_routed.vcf.gz",
        index="results/variants/stage1_routed.vcf.gz.csi"

    output:
        bcf="results/variants/stage2_variant_keys.bcf",
        index="results/variants/stage2_variant_keys.bcf.csi"

    shell:
        r"""
        set -euo pipefail

        TMP="{output.bcf}.tmp.bcf"
        MAP="{output.bcf}.chrmap.tmp"

        rm -f "$TMP" "$TMP.csi" "$MAP"

        for CHR in $(seq 1 22) X Y; do
            printf '%s\tchr%s\n' "$CHR" "$CHR"
        done > "$MAP"

        printf 'MT\tchrM\n' >> "$MAP"
        printf 'M\tchrM\n' >> "$MAP"

        bcftools view \
            -G \
            {input.vcf} \
            -Ou \
        | bcftools annotate \
            -x ID,FILTER,INFO \
            --rename-chrs "$MAP" \
            -Ob \
            -o "$TMP"

        bcftools index -f "$TMP"

        test "$(bcftools index -n "$TMP")" -gt 0

        mv "$TMP" {output.bcf}
        mv "$TMP.csi" {output.index}

        rm -f "$MAP"
        """

GNOMAD_NUCLEAR_CHROMS = [
    str(i) for i in range(1, 23)
] + ["X", "Y"]

GNOMAD_JOINT_URL = config["stage2"]["gnomad"]["url_template"]

GNOMAD_MT = config["stage2"]["gnomad"]["mt_vcf"]

rule stage2_gnomad_remote_chr:
    conda:
        "../../envs/environment.yaml"
    input:
        variant_keys="results/variants/stage2_variant_keys.bcf",
        variant_keys_index="results/variants/stage2_variant_keys.bcf.csi"

    output:
        "results/intermediate/stage2/gnomad_remote_by_chr/chr{chrom}.parquet"

    params:
        url=lambda wc: GNOMAD_JOINT_URL.format(chrom=wc.chrom),
        retries=4,
        retry_delay_seconds=30,
        remote_timeout_seconds=config["stage2"].get("gnomad_remote_timeout_seconds", 7200),
        index_cache_dir=config["stage2"]["gnomad"].get(
            "index_cache_dir", "data/resources/gnomad/joint/v4.1/indexes"
        )

    wildcard_constraints:
        chrom="(?:[1-9]|1[0-9]|2[0-2]|X|Y)"

    threads: 1

    resources:
        gnomad_remote=1

    script:
        "../scripts/extract_stage2_gnomad_remote_chr.py"

    
rule stage2_add_gnomad:
    conda:
        "../../envs/environment.yaml"
    input:
        stage2="results/tables/stage2_safety_complete.parquet",

        nuclear=expand(
            "results/intermediate/stage2/gnomad_remote_by_chr/chr{chrom}.parquet",
            chrom=GNOMAD_NUCLEAR_CHROMS,
        ),

        mt_vcf=GNOMAD_MT,
        mt_index=GNOMAD_MT + ".tbi"

    output:
        "results/tables/stage2_safety_annotated.parquet"

    script:
        "../scripts/merge_stage2_gnomad_remote_v2.py"
