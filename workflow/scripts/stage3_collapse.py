from pathlib import Path
import os

import pandas as pd
from workflow_checks import validate_table, check_transition


# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

SEP = " | "

KEY = [
    "Chr_Pos",
    "Variant_REF_ALT",
]

# These fields are transcript-aligned.
# Their order must remain identical after collapse.
TX_FIELDS = [
    "Transcript_ID",
    "Gene_symbol",
    "Gene_ID",
    "Transcript_biotype",
    "Consequence",
    "Genomic_region",
]

# These should already be identical for every transcript row
# belonging to the same variant.
VARIANT_FIELDS = [
    "FILTER", "GT", "GQ", "DP", "AD_ref", "AD_alt", "VAF",

    # Nuclear gnomAD
    "gnomAD_AF",
    "gnomAD_grpmax_AF",
    "gnomAD_grpmax_population",
    "gnomAD_het",
    "gnomAD_hom",
    "gnomAD_hemi",

    # Mitochondrial gnomAD
    "gnomAD_mt_AN",
    "gnomAD_mt_AC_hom",
    "gnomAD_mt_AC_het",
    "gnomAD_mt_AF_hom",
    "gnomAD_mt_AF_het",
    "gnomAD_mt_max_hl",

    # ClinVar
    "ClinVar_ID",
    "ClinVar_significance",
]


# VEP consequence severity order.
# Used only to derive Most_severe_consequence.
CONSEQUENCE_ORDER = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "feature_elongation",
    "feature_truncation",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_donor_5th_base_variant",
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "regulatory_region_variant",
    "intergenic_variant",
]

CONSEQUENCE_RANK = {
    consequence: rank
    for rank, consequence in enumerate(CONSEQUENCE_ORDER)
}


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def missing(x):
    return pd.isna(x) or str(x).strip() in {
        "",
        ".",
        "nan",
        "None",
        "<NA>",
    }


def clean(x):
    """
    Preserve transcript alignment.

    Missing transcript-level values become '-'
    instead of disappearing from the joined string.
    """
    return "-" if missing(x) else str(x).strip()


def yes(x):
    """
    Interpret MANE_Select-like fields.
    """
    return (
        not missing(x)
        and str(x).strip().lower()
        not in {"0", "false", "no"}
    )


def join_unique(values):
    """
    Join unique non-missing values while preserving order.
    """
    return SEP.join(
        dict.fromkeys(
            str(x).strip()
            for x in values
            if not missing(x)
        )
    )


def most_severe(values):
    """
    Find the most severe VEP consequence across all transcripts.
    """
    terms = []

    for value in values:
        if missing(value):
            continue

        terms.extend(
            term.strip()
            for term in str(value).split("&")
            if term.strip()
        )

    return min(
        terms,
        key=lambda term: CONSEQUENCE_RANK.get(term, 10**6),
        default="",
    )


# ---------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------

def load_stage2(path):
    required_columns = (
        KEY
        + TX_FIELDS
        + VARIANT_FIELDS
        + [
            "MANE_Select",
            "Regulatory_context",
        ]
    )

    try:
        return pd.read_parquet(
            path,
            columns=required_columns,
        )

    except Exception as exc:
        raise ValueError(
            "Could not read the required Stage 2 columns. "
            "Check the Stage 2 parquet schema."
        ) from exc


# ---------------------------------------------------------------------
# Collapse variant-level scalar fields
# ---------------------------------------------------------------------

def collapse_variant_fields(df):
    """
    Collapse fields that should be identical across transcript rows.

    One safety check is retained:
    if a supposedly variant-level value differs between transcripts
    of the same variant, stop rather than silently selecting one.
    """

    grouped = df.groupby(KEY, sort=False)

    conflicts = (
        grouped[VARIANT_FIELDS]
        .nunique(dropna=True)
        .gt(1)
        .any(axis=1)
    )

    if conflicts.any():
        examples = conflicts[conflicts].index.tolist()[:5]

        raise ValueError(
            "Variant-level values disagree across transcript rows "
            f"for {int(conflicts.sum()):,} variants.\n"
            f"Examples: {examples}"
        )

    return (
        grouped[VARIANT_FIELDS]
        .first()
        .reset_index()
    )


# ---------------------------------------------------------------------
# Collapse transcript annotations
# ---------------------------------------------------------------------

def collapse_transcripts(df):
    """
    Collapse transcript annotations while keeping fields aligned.

    Transcript ordering:
        1. MANE Select
        2. remaining transcripts sorted by Transcript_ID
    """

    tx = df[
        KEY
        + TX_FIELDS
        + ["MANE_Select"]
    ].copy()

    tx["_mane"] = tx["MANE_Select"].map(yes)

    tx["_priority"] = 1
    tx.loc[tx["_mane"], "_priority"] = 0

    tx["_tx_sort"] = tx["Transcript_ID"].map(clean)

    tx = tx.sort_values(
        KEY + ["_priority", "_tx_sort"],
        kind="stable",
    )

    # Same transcript may occasionally appear more than once.
    tx = tx.drop_duplicates(
        KEY + ["Transcript_ID"],
        keep="first",
    )

    # Transcript-aligned joined fields.
    agg = {
        field: (
            lambda values:
            SEP.join(clean(x) for x in values)
        )
        for field in TX_FIELDS
    }

    collapsed = (
        tx.groupby(KEY, sort=False)
        .agg(agg)
        .reset_index()
    )

    collapsed = collapsed.rename(
        columns={
            "Consequence": "Consequences",
        }
    )

    # MANE transcript IDs.
    mane = (
        tx[tx["_mane"]]
        .groupby(KEY)["Transcript_ID"]
        .agg(join_unique)
        .reset_index(
            name="MANE_Transcript_ID"
        )
    )

    # Most severe consequence.
    severe = (
        tx.groupby(KEY)["Consequence"]
        .agg(most_severe)
        .reset_index(
            name="Most_severe_consequence"
        )
    )

    collapsed = (
        collapsed
        .merge(
            mane,
            on=KEY,
            how="left",
        )
        .merge(
            severe,
            on=KEY,
            how="left",
        )
    )

    collapsed["MANE_Transcript_ID"] = (
        collapsed["MANE_Transcript_ID"]
        .fillna("")
    )

    return collapsed


# ---------------------------------------------------------------------
# Regulatory annotation
# ---------------------------------------------------------------------

def collapse_regulatory(df):
    """
    Regulatory context is variant-level information.

    Keep all unique regulatory annotations.
    It is intentionally NOT transcript-aligned.
    """

    return (
        df.groupby(KEY, sort=False)["Regulatory_context"]
        .agg(join_unique)
        .reset_index()
    )


# ---------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------

def validate(final):
    """
    Confirm:
      - one row per variant
      - transcript fields remain aligned
    """

    if final[KEY].duplicated().any():
        raise ValueError(
            "Duplicate variants remain after Stage 3 collapse."
        )

    aligned_fields = [
        "Transcript_ID",
        "Gene_symbol",
        "Gene_ID",
        "Transcript_biotype",
        "Consequences",
        "Genomic_region",
    ]

    def aligned(row):
        lengths = {
            len(str(row[field]).split(SEP))
            for field in aligned_fields
        }
        return len(lengths) == 1

    bad = ~final.apply(
        aligned,
        axis=1,
    )

    if bad.any():
        example = final.loc[bad].iloc[0]

        raise ValueError(
            "Transcript fields lost alignment at "
            f"{example['Chr_Pos']} "
            f"{example['Variant_REF_ALT']}"
        )


# ---------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------

def save_parquet(df, path):
    """
    Atomic output:
    write temporary file first, then replace final path.
    """

    path = Path(path)
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    tmp = Path(str(path) + ".tmp")

    df.to_parquet(
        tmp,
        index=False,
        compression="snappy",
    )

    os.replace(tmp, path)


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    stage2_path = str(
        snakemake.input.stage2
    )

    output_path = str(
        snakemake.output.collapsed
    )

    print(f"Reading Stage 2: {stage2_path}")

    df = load_stage2(stage2_path)

    print(
        f"Stage 2 rows: "
        f"{len(df):,}"
    )

    print(
        "Collapsing variant-level fields..."
    )

    validate_table(df, "Stage 2 to Stage 3", unique=False)
    repeated = df.duplicated(KEY + ["Transcript_ID"], keep=False)
    if repeated.any():
        groups = df.loc[repeated].groupby(KEY + ["Transcript_ID"], dropna=False)
        if groups[TX_FIELDS[1:]].nunique(dropna=True).gt(1).any().any():
            raise ValueError("Conflicting gene/transcript annotations would be lost during collapse")

    variants = collapse_variant_fields(df)

    print(
        "Collapsing transcript annotations..."
    )

    transcripts = collapse_transcripts(df)

    print(
        "Collapsing regulatory annotations..."
    )

    regulatory = collapse_regulatory(df)

    print(
        "Building final collapsed table..."
    )

    final = (
        variants
        .merge(
            transcripts,
            on=KEY,
            validate="one_to_one",
        )
        .merge(
            regulatory,
            on=KEY,
            how="left",
            validate="one_to_one",
        )
    )

    final["Regulatory_context"] = (
        final["Regulatory_context"]
        .fillna("")
    )

    validate(final)
    check_transition(df, final, "stage3_collapse", operation="collapse", unique=False, aligned=True, required=["MANE_Transcript_ID"])

    save_parquet(
        final,
        output_path,
    )

    print("\nStage 3 collapse complete")
    print("-" * 45)
    print(
        f"Input transcript rows: "
        f"{len(df):,}"
    )
    print(
        f"Unique variants:       "
        f"{len(final):,}"
    )
    print(
        f"Output: {output_path}"
    )


if __name__ == "__main__":
    main()