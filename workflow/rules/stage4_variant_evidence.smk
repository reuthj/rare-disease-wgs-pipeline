rule stage4_route:
    conda:
        "../../envs/environment.yaml"
    input:
        candidates="results/tables/stage3_candidates.csv"

    output:
        candidates="results/intermediate/stage4/stage4_variants_routed.csv"

    script:
        "../scripts/stage4_route.py"


rule stage4_variant_evidence:
    input:
        variants="results/intermediate/stage4/stage4_variants_routed.csv",
        vep_cache="results/intermediate/stage3/vep_after_rarity.json.gz",
        gtf=config["resources"]["gencode"],
        script="workflow/scripts/stage4_add_evidence.py",
    output:
        "results/tables/stage4_variant_evidence.csv",
    params:
        spliceai_conda_env=config["stage4"].get(
            "spliceai_conda_env",
            "spliceai-env",
        ),
        spliceai_conda_prefix=config["stage4"].get("spliceai_conda_prefix") or "",
    conda:
        "../../envs/environment.yaml"
    shell:
        """
        python {input.script} \
            --input {input.variants} \
            --output {output} \
            --vep-cache {input.vep_cache} \
            --gtf {input.gtf} \
            --stage4-dir results/intermediate/stage4 \
            --clinvar-dir data/resources/clinvar \
            --alphamissense-dir data/resources/alphamissense \
            --spliceai-conda-env {params.spliceai_conda_env:q} \
            --spliceai-conda-prefix={params.spliceai_conda_prefix:q}
        """
