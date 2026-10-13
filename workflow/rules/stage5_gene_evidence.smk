rule stage5_gene_evidence:
    conda:
        "../../envs/environment.yaml"
    input:
        variants="results/tables/stage4_variant_evidence.csv",
        case_context="config/case_context.yaml",
        gencc=config["stage5"]["resources"]["gencc_submissions"],
        constraint=config["stage5"]["resources"]["gnomad_constraint"],
        clinvar=config["stage5"]["resources"]["clinvar_vcf"],
        clinvar_index=config["stage5"]["resources"]["clinvar_index"],
        hpo_annotations=config["stage5"]["resources"]["hpo_annotations"],
        hpo_ontology=config["stage5"]["resources"]["hpo_ontology"],
        runner="workflow/scripts/stage5_gene_evidence.py",
    output:
        variants="results/intermediate/stage5/stage5_gene_evidence.csv",
    params:
        cache_dir="results/intermediate/stage5/cache",
        progress_every=500000,
    shell:
        """
        python {input.runner:q} \
            --input {input.variants:q} \
            --output {output.variants:q} \
            --case-context {input.case_context:q} \
            --gencc {input.gencc:q} \
            --constraint {input.constraint:q} \
            --clinvar {input.clinvar:q} \
            --hpo-annotations {input.hpo_annotations:q} \
            --hpo-ontology {input.hpo_ontology:q} \
            --cache-dir {params.cache_dir:q} \
            --progress-every {params.progress_every}
        """
        
rule stage5_literature:
    conda:
        "../../envs/environment.yaml"
    input:
        variants="results/intermediate/stage5/stage5_gene_evidence.csv",
        case_context="config/case_context.yaml",
        script="workflow/scripts/stage5_literature.py",

    output:
        variants="results/tables/stage5_gene_evidence_complete.csv",

    params:
        cache_dir=config["stage5"]["literature"].get(
            "cache_dir",
            "results/intermediate/stage5/cache/europe_pmc",
        ),
        progress_every=config["stage5"]["literature"].get(
            "progress_every",
            25,
        ),
        max_retries=config["stage5"]["europe_pmc"].get("max_retries", 3),

    shell:
        """
        python {input.script:q} \
            --input {input.variants:q} \
            --output {output.variants:q} \
            --case-context {input.case_context:q} \
            --cache-dir {params.cache_dir:q} \
            --progress-every {params.progress_every} \
            --max-retries {params.max_retries}
        """

    
