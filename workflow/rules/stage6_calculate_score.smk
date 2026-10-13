rule stage6_score:
    input:
        variants=config["stage6"]["input_table"],
        case_context=config["stage6"]["case_context"],
    output:
        variants=config["stage6"]["output_table"],
    params:
        progress_every=config["stage6"].get("progress_every", 100),
    conda:
        "../../envs/environment.yaml"
    script:
        "../scripts/stage6_score.py"