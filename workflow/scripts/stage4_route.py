import re
from pathlib import Path

import pandas as pd
from workflow_checks import validate_table, check_transition


# ============================================================
# Settings
# ============================================================

TX_SEP = " | "
MISSING = {"", ".", "-", "_", "nan", "none", "<na>"}

SPLICE_MAX = 50

MISSENSE = {
    "missense_variant",
}

LOF = {
    "frameshift_variant",
    "stop_gained",
    "splice_acceptor_variant",
    "splice_donor_variant",
}

INFRAME = {
    "inframe_insertion",
    "inframe_deletion",
}

SPLICE = {
    "splice_acceptor_variant",
    "splice_donor_variant",
    "splice_region_variant",
    "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
}


# ============================================================
# Small helpers
# ============================================================

def is_missing(value):
    if pd.isna(value):
        return True
    return str(value).strip().lower() in MISSING


def strip_version(value):
    if is_missing(value):
        return ""
    return str(value).strip().split(".")[0]


def split_field(value):
    if is_missing(value):
        return []
    return [
        x.strip()
        for x in str(value).split(TX_SEP)
    ]


def consequence_terms(value):
    if is_missing(value):
        return set()

    return {
        x.strip()
        for x in re.split(r"[|&,]", str(value))
        if x.strip()
    }


def parse_tx_map(value):
    """
    ENST1=exon 3/10 | ENST2=intron 2/8
    ->
    {
        ENST1: exon 3/10,
        ENST2: intron 2/8
    }
    """

    result = {}

    for item in split_field(value):

        if "=" not in item:
            continue

        tx, info = item.split("=", 1)

        tx = strip_version(tx)

        if tx:
            result[tx] = info.strip()

    return result


# ============================================================
# Main
# ============================================================

df = pd.read_csv(
    snakemake.input.candidates,
    low_memory=False,
)

print("=" * 60)
print("STAGE 4.1 — VARIANT ROUTING")
print("=" * 60)
print(f"Input variants: {len(df):,}")


validate_table(df, "Stage 4 route input", required=["MANE_Transcript_ID"], aligned=True)
route_rows = []

for _, row in df.iterrows():

    txs = split_field(
        row["Transcript_ID"]
    )

    consequences = split_field(
        row["Consequences"]
    )

    if len(txs) != len(consequences):
        raise ValueError(
            f"Transcript/consequence alignment failed: "
            f"{row['Chr_Pos']} "
            f"{row['Variant_REF_ALT']}"
        )

    gene_context = parse_tx_map(
        row["Gene_context"]
    )

    splice_context = parse_tx_map(
        row["Splice_context"]
    )

    tx_routes = {}
    variant_routes = set()

    # --------------------------------------------------------
    # Transcript-level routing
    # --------------------------------------------------------

    for tx, consequence in zip(
        txs,
        consequences,
    ):

        tx = strip_version(tx)

        if not tx:
            continue

        csq = consequence_terms(
            consequence
        )

        context = (
            gene_context
            .get(tx, "")
            .lower()
        )

        splice_distance = pd.to_numeric(
            splice_context.get(tx),
            errors="coerce",
        )

        routes = set()

        if (
            context.startswith("exon")
            and csq & MISSENSE
        ):
            routes.add("MISSENSE")

        if csq & LOF:
            routes.add("LOF")

        if (
            context.startswith("exon")
            and csq & INFRAME
        ):
            routes.add("INFRAME_INDEL")

        if csq & SPLICE:
            routes.add("SPLICE")

        elif (
            context.startswith(
                ("exon", "intron")
            )
            and pd.notna(
                splice_distance
            )
            and abs(
                float(splice_distance)
            ) <= SPLICE_MAX
        ):
            routes.add("SPLICE")

        tx_routes[tx] = routes
        variant_routes.update(routes)

    # --------------------------------------------------------
    # Variant-level routes already established by Stage 3
    # --------------------------------------------------------

    if not is_missing(
        row.get("Regulatory_context")
    ):
        variant_routes.add(
            "REGULATORY"
        )

    chrom = str(
        row["Chr_Pos"]
    ).split(":")[0].removeprefix("chr")

    if chrom in {"M", "MT"}:
        variant_routes.add(
            "MITOCHONDRIAL"
        )

    # --------------------------------------------------------
    # Only expose transcript routes when transcripts differ
    # --------------------------------------------------------

    transcript_routes = pd.NA

    route_profiles = {
        frozenset(routes)
        for routes in tx_routes.values()
    }

    if len(route_profiles) > 1:

        transcript_routes = TX_SEP.join(
            f"{tx}:"
            + (
                ",".join(sorted(routes))
                if routes
                else "COMMON"
            )
            for tx, routes
            in tx_routes.items()
        )

    route_rows.append({

        "is_missense":
            "MISSENSE"
            in variant_routes,

        "is_lof":
            "LOF"
            in variant_routes,

        "is_inframe_indel":
            "INFRAME_INDEL"
            in variant_routes,

        "is_splice_relevant":
            "SPLICE"
            in variant_routes,

        "is_regulatory":
            "REGULATORY"
            in variant_routes,

        "is_mitochondrial":
            "MITOCHONDRIAL"
            in variant_routes,

        "Evidence_routes":
            (
                "COMMON"
                if not variant_routes
                else
                "COMMON | "
                + " | ".join(
                    sorted(
                        variant_routes
                    )
                )
            ),

        "Transcript_routes":
            transcript_routes,
    })


routes = pd.DataFrame(
    route_rows,
    index=df.index,
)

overlap = (
    set(df.columns)
    & set(routes.columns)
)

if overlap:
    raise RuntimeError(
        "Stage 4.1 attempted to overwrite "
        "existing columns: "
        + ", ".join(sorted(overlap))
    )

final = pd.concat(
    [df, routes],
    axis=1,
)


# ============================================================
# Validation
# ============================================================

if len(final) != len(df):
    raise RuntimeError(
        "Stage 4.1 changed variant count."
    )


print("\nEvidence routes")
print("-" * 60)

for col in [
    "is_missense",
    "is_lof",
    "is_inframe_indel",
    "is_splice_relevant",
    "is_regulatory",
    "is_mitochondrial",
]:
    print(
        f"{col:<25} "
        f"{int(final[col].sum()):>8,}"
    )


Path(
    snakemake.output.candidates
).parent.mkdir(
    parents=True,
    exist_ok=True,
)

check_transition(df, final, "stage4_route", aligned=True)

final.to_csv(
    snakemake.output.candidates,
    index=False,
)

print("\nValidation")
print("-" * 60)
print(f"Input variants:  {len(df):,}")
print(f"Output variants: {len(final):,}")
print("Stage 4.1: PASS")