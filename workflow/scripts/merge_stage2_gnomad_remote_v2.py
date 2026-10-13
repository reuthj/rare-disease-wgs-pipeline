import os
import subprocess
from io import StringIO

import pandas as pd
from workflow_checks import validate_table, check_transition


INPUT_PARQUET = str(snakemake.input.stage2)
NUCLEAR_FILES = [str(path) for path in snakemake.input.nuclear]
MT_VCF = str(snakemake.input.mt_vcf)
OUTPUT_PARQUET = str(snakemake.output[0])

CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]
KEY = ["_chrom", "_pos", "_ref", "_alt"]

# Only the compact, clinically useful nuclear gnomAD fields enter final Stage 2.
NUCLEAR_FINAL_FIELDS = [
    "gnomAD_AF",
    "gnomAD_grpmax_AF",
    "gnomAD_grpmax_population",
    "gnomAD_AC",
    "gnomAD_AN",
    "gnomAD_het",
    "gnomAD_hom",
    "gnomAD_hemi",
]

# gnomAD v3.1 mtDNA fields kept in final Stage 2.
MT_FIELDS = {
    "AN": ("gnomAD_mt_AN", "Int64"),
    "AC_hom": ("gnomAD_mt_AC_hom", "Int64"),
    "AC_het": ("gnomAD_mt_AC_het", "Int64"),
    "AF_hom": ("gnomAD_mt_AF_hom", "float"),
    "AF_het": ("gnomAD_mt_AF_het", "float"),
    "max_hl": ("gnomAD_mt_max_hl", "float"),
}
MT_FINAL_FIELDS = [output for output, _ in MT_FIELDS.values()]


def canon_chrom(chrom):
    chrom = str(chrom).removeprefix("chr")
    return "MT" if chrom in {"M", "MT"} else chrom


def add_key(df):
    chrom_pos = df["Chr_Pos"].astype(str).str.rsplit(":", n=1, expand=True)
    ref_alt = df["Variant_REF_ALT"].astype(str).str.split(">", n=1, expand=True)

    out = df.copy()
    out["_chrom"] = chrom_pos[0].map(canon_chrom)
    out["_pos"] = pd.to_numeric(chrom_pos[1], errors="raise").astype("int64")
    out["_ref"] = ref_alt[0]
    out["_alt"] = ref_alt[1]
    return out


def load_nuclear_checkpoints():
    frames = []

    for path in NUCLEAR_FILES:
        df = pd.read_parquet(path)

        required = set(KEY + NUCLEAR_FINAL_FIELDS)
        missing = required.difference(df.columns)
        if missing:
            raise ValueError(
                f"Nuclear checkpoint is missing required columns: {path}: "
                f"{sorted(missing)}"
            )

        # Keep only fields owned by final Stage 2. Raw/helper gnomAD fields,
        # if present in a checkpoint, are intentionally not propagated.
        frames.append(df[KEY + NUCLEAR_FINAL_FIELDS].copy())

    nuclear = pd.concat(frames, ignore_index=True)

    if nuclear.duplicated(KEY).any():
        raise ValueError("Duplicate variant keys across nuclear gnomAD checkpoints")

    nuclear["_checkpoint_seen"] = True
    return nuclear


def check_mt_header():
    header = subprocess.run(
        ["bcftools", "view", "-h", MT_VCF],
        check=True,
        capture_output=True,
        text=True,
    ).stdout

    missing = [tag for tag in MT_FIELDS if f"##INFO=<ID={tag}," not in header]
    if missing:
        raise ValueError(f"Missing mtDNA INFO fields: {missing}")


def query_mt(targets):
    """Read the tiny local gnomAD v3.1 chrM VCF once and keep project mtDNA matches."""
    mt_targets = targets.loc[targets["_chrom"] == "MT", KEY].drop_duplicates()
    if mt_targets.empty:
        return pd.DataFrame(columns=KEY + ["_mt_found"] + MT_FINAL_FIELDS)

    check_mt_header()

    tags = list(MT_FIELDS)
    fmt = (
        "%CHROM\t%POS\t%REF\t%ALT\t"
        + "\t".join(f"%INFO/{tag}" for tag in tags)
        + "\n"
    )
    names = ["CHROM", "POS", "REF", "ALT", *tags]

    result = subprocess.run(
        ["bcftools", "query", "-f", fmt, MT_VCF],
        check=True,
        capture_output=True,
        text=True,
    )

    if not result.stdout.strip():
        return pd.DataFrame(columns=KEY + ["_mt_found"] + MT_FINAL_FIELDS)

    g = pd.read_csv(
        StringIO(result.stdout),
        sep="\t",
        header=None,
        names=names,
        na_values=".",
        keep_default_na=True,
        dtype="string",
    )

    g["_chrom"] = "MT"
    g["_pos"] = pd.to_numeric(g["POS"], errors="raise").astype("int64")
    g["_ref"] = g["REF"]
    g["_alt"] = g["ALT"]

    # Exact variant match against the Stage-2 mitochondrial target set.
    g = g.merge(mt_targets, on=KEY, how="inner", validate="many_to_one")

    g = g.rename(
        columns={tag: output for tag, (output, _) in MT_FIELDS.items()}
    )

    for _, (column, dtype) in MT_FIELDS.items():
        if dtype == "Int64":
            g[column] = pd.to_numeric(g[column], errors="coerce").astype("Int64")
        elif dtype == "float":
            g[column] = pd.to_numeric(g[column], errors="coerce")

    g["_mt_found"] = True

    if g.duplicated(KEY).any():
        raise ValueError("Duplicate mitochondrial gnomAD variant keys")

    return g[KEY + ["_mt_found"] + MT_FINAL_FIELDS]


def main():
    stage2 = add_key(pd.read_parquet(INPUT_PARQUET))
    original = stage2.drop(columns=KEY).copy()
    validate_table(original, "Stage 2 population input", unique=False)
    original_rows = len(stage2)
    variants = stage2[KEY].drop_duplicates()

    nuclear = load_nuclear_checkpoints()

    # Clean rerun: replace only the compact fields owned by this stage.
    stage2 = stage2.drop(
        columns=[
            column
            for column in NUCLEAR_FINAL_FIELDS + MT_FINAL_FIELDS
            if column in stage2.columns
        ],
        errors="ignore",
    )

    out = stage2.merge(
        nuclear,
        on=KEY,
        how="left",
        validate="many_to_one",
    )

    if len(out) != original_rows:
        raise ValueError(
            f"Row count changed during nuclear merge: {original_rows} -> {len(out)}"
        )

    # Every nuclear Stage-2 variant must appear in exactly one completed checkpoint.
    nuclear_rows = out["_chrom"].isin(CHROMS)
    missing_checkpoint = nuclear_rows & out["_checkpoint_seen"].isna()

    if missing_checkpoint.any():
        examples = (
            out.loc[missing_checkpoint, ["Chr_Pos", "Variant_REF_ALT"]]
            .drop_duplicates()
            .head(10)
        )
        missing_n = out.loc[missing_checkpoint, KEY].drop_duplicates().shape[0]
        raise ValueError(
            f"{missing_n:,} nuclear variants are missing from chromosome checkpoints. "
            "Examples:\n" + examples.to_string(index=False)
        )

    # mtDNA is tiny, so query it locally here instead of creating another job/checkpoint.
    mt = query_mt(variants)
    out = out.merge(
        mt,
        on=KEY,
        how="left",
        validate="many_to_one",
    )

    if len(out) != original_rows:
        raise ValueError(
            f"Row count changed during mtDNA merge: {original_rows} -> {len(out)}"
        )

    # A project mtDNA variant absent from the gnomAD mt sites VCF has zero observed
    # homoplasmic/heteroplasmic alternate counts. AN/AF/max_hl stay missing because
    # the absent site has no gnomAD row from which to obtain those values.
    mt_rows = out["_chrom"] == "MT"
    mt_absent = mt_rows & out["_mt_found"].isna()

    for column in ["gnomAD_mt_AC_hom", "gnomAD_mt_AC_het"]:
        out[column] = pd.to_numeric(out[column], errors="coerce").astype("Int64")
        out.loc[mt_absent, column] = 0

    out = out.drop(columns=KEY + ["_checkpoint_seen", "_mt_found"])

    os.makedirs(os.path.dirname(OUTPUT_PARQUET), exist_ok=True)
    tmp = OUTPUT_PARQUET + ".tmp"
    check_transition(original, out, "stage2_population", unique=False, changed_fields=NUCLEAR_FINAL_FIELDS + MT_FINAL_FIELDS)
    out.to_parquet(tmp, index=False)
    os.replace(tmp, OUTPUT_PARQUET)

    nuclear_unique = variants[variants["_chrom"].isin(CHROMS)].shape[0]
    mt_unique = variants[variants["_chrom"] == "MT"].shape[0]

    print(f"Stage-2 rows preserved: {len(out):,}")
    print(f"Unique nuclear variants: {nuclear_unique:,}")
    print(f"Unique mtDNA variants: {mt_unique:,}")
    print(f"Nuclear gnomAD fields kept: {len(NUCLEAR_FINAL_FIELDS)}")
    print(f"mtDNA gnomAD fields kept: {len(MT_FINAL_FIELDS)}")
    print(f"Output: {OUTPUT_PARQUET}")


if __name__ == "__main__":
    main()
