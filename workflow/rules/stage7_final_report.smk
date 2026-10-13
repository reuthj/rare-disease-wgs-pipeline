rule stage7_final_report:
    conda:
        "../../envs/environment.yaml"
    input:
        variants=config["stage7"]["input_table"],
        case_context=config["stage7"]["case_context"],
    output:
        workbook=config["stage7"]["output_workbook"],
    script:
        "../scripts/stage7_final_report.py"