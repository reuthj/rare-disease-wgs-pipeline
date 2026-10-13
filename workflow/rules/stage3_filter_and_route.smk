rule stage3_collapse:
    conda:
        "../../envs/environment.yaml"
    input:
        stage2="results/tables/stage2_safety_annotated.parquet"
    output:
        collapsed="results/intermediate/stage3/stage3_collapsed.parquet"
    script:
        "../scripts/stage3_collapse.py"

rule stage3_filter_route_reviewed:
    conda:
        "../../envs/environment.yaml"
    input:
        collapsed="results/intermediate/stage3/stage3_collapsed.parquet",
        case_context="config/case_context.yaml",
        acmg_sf=config["resources"]["acmg_sf"]["file"],
        gencode_gtf=config["resources"]["gencode"],
    output:
        candidates="results/tables/stage3_candidates.csv",
        filtered_keys="results/variants/stage3_filtered_keys.tsv.gz",
        vep_cache="results/intermediate/stage3/vep_after_rarity.json.gz",
    script:
        "../scripts/stage3_filter_route.py"
