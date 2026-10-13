import os
import subprocess
import tempfile

import pandas as pd
from workflow_checks import validate_table, check_transition


INPUT_PARQUET = snakemake.input.parquet
INPUT_VCF = snakemake.input.vcf
OUTPUT_PARQUET = snakemake.output[0]


def canon_chrom(chrom):
    chrom = str(chrom)
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    return "MT" if chrom == "M" else chrom


def split_variant_fields(df):
    chrom_pos = df["Chr_Pos"].astype(str).str.rsplit(":", n=1, expand=True)
    ref_alt = df["Variant_REF_ALT"].astype(str).str.split(">", n=1, expand=True)

    out = df.copy()
    out["_chrom"] = chrom_pos[0].map(canon_chrom)
    out["_pos"] = pd.to_numeric(chrom_pos[1], errors="raise").astype("int64")
    out["_ref"] = ref_alt[0]
    out["_alt"] = ref_alt[1]
    return out


def vcf_contig_map(vcf):
    result = subprocess.run(
        ["bcftools", "index", "-s", vcf],
        check=True,
        capture_output=True,
        text=True,
    )
    contigs = [line.split("\t", 1)[0] for line in result.stdout.splitlines() if line]
    return {canon_chrom(contig): contig for contig in contigs}


def parse_ad(value):
    if pd.isna(value) or value in {"", "."}:
        return pd.NA, pd.NA
    parts = str(value).split(",")
    if len(parts) < 2:
        return pd.NA, pd.NA
    ref = pd.to_numeric(parts[0], errors="coerce")
    alt = pd.to_numeric(parts[1], errors="coerce")
    return ref, alt


def main():
    df = pd.read_parquet(INPUT_PARQUET)
    validate_table(df, "Stage 2 genotype input", unique=False)
    original = df
    original_rows = len(df)
    df = split_variant_fields(df)

    variants = df[["_chrom", "_pos", "_ref", "_alt"]].drop_duplicates().copy()
    contig_map = vcf_contig_map(INPUT_VCF)
    variants["_vcf_chrom"] = variants["_chrom"].map(contig_map)

    missing_contigs = variants.loc[variants["_vcf_chrom"].isna(), "_chrom"].unique()
    if len(missing_contigs):
        raise ValueError(f"Stage 2 chromosomes not found in normalized VCF: {missing_contigs.tolist()}")

    samples = subprocess.run(
        ["bcftools", "query", "-l", INPUT_VCF],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.splitlines()
    if len(samples) != 1:
        raise ValueError(f"Expected one proband sample in normalized VCF, found {len(samples)}")

    with tempfile.TemporaryDirectory() as tmpdir:
        regions_path = os.path.join(tmpdir, "stage2.regions")
        query_path = os.path.join(tmpdir, "genotypes.tsv")

        regions = variants[["_vcf_chrom", "_pos"]].drop_duplicates().sort_values(["_vcf_chrom", "_pos"])
        with open(regions_path, "w") as handle:
            for chrom, pos in regions.itertuples(index=False):
                handle.write(f"{chrom}\t{pos}\t{pos}\n")

        with open(query_path, "w") as handle:
            subprocess.run(
                [
                    "bcftools", "query",
                    "-R", regions_path,
                    "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER[\\t%GT\\t%DP\\t%GQ\\t%AD]\\n",
                    INPUT_VCF,
                ],
                check=True,
                stdout=handle,
                text=True,
            )

        geno = pd.read_csv(
            query_path,
            sep="\t",
            header=None,
            names=["_vcf_chrom","_pos","_ref","_alt","FILTER","GT","DP","GQ","AD",],
            dtype={"_vcf_chrom": "string","_ref": "string","_alt": "string","FILTER": "string","GT": "string","AD": "string",
            },
        )

    geno["_chrom"] = geno["_vcf_chrom"].map(canon_chrom)
    geno["_pos"] = pd.to_numeric(geno["_pos"], errors="raise").astype("int64")

    key = ["_chrom", "_pos", "_ref", "_alt"]
    if geno.duplicated(key).any():
        raise ValueError("Normalized VCF returned duplicate records for a Stage 2 variant key")

    ad = geno["AD"].apply(parse_ad)
    geno["AD_ref"] = [x[0] for x in ad]
    geno["AD_alt"] = [x[1] for x in ad]
    geno["DP"] = pd.to_numeric(geno["DP"], errors="coerce").astype("Int64")
    geno["GQ"] = pd.to_numeric(geno["GQ"], errors="coerce").astype("Int64")
    geno["AD_ref"] = pd.to_numeric(geno["AD_ref"], errors="coerce").astype("Int64")
    geno["AD_alt"] = pd.to_numeric(geno["AD_alt"], errors="coerce").astype("Int64")


    denom = geno["AD_ref"] + geno["AD_alt"]
    geno["VAF"] = (geno["AD_alt"] / denom).where(denom > 0)
    geno["_matched"] = True

    keep = key + ["FILTER","GT", "DP","GQ", "AD_ref", "AD_alt", "VAF", "_matched"]
    out = df.merge(geno[keep], on=key, how="left", validate="many_to_one")

    if len(out) != original_rows:
        raise ValueError(f"Row count changed during genotype merge: {original_rows} -> {len(out)}")

    unmatched = out.loc[out["_matched"].isna(), ["Chr_Pos", "Variant_REF_ALT"]].drop_duplicates()
    if len(unmatched):
        examples = unmatched.head(10).to_dict("records")
        raise ValueError(f"{len(unmatched)} Stage 2 variants were not found in normalized VCF. Examples: {examples}")

    out = out.drop(columns=key + ["_matched"])

    os.makedirs(os.path.dirname(OUTPUT_PARQUET), exist_ok=True)
    tmp_output = OUTPUT_PARQUET + ".tmp"
    check_transition(original, out, "stage2_genotype", unique=False)
    out.to_parquet(tmp_output, index=False)
    os.replace(tmp_output, OUTPUT_PARQUET)

    print(f"Stage 2 rows: {original_rows:,}")
    print(f"Unique variants: {len(variants):,}")
    print(f"Output: {OUTPUT_PARQUET}")
    print("Added: FILTER, GT, DP, GQ, AD_ref, AD_alt, VAF")


if __name__ == "__main__":
    main()
