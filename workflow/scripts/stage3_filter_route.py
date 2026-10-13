#!/usr/bin/env python3

from collections import defaultdict
import gzip
import re
import time
import json
from pathlib import Path
import requests
from main_classification import is_main_plp
import pandas as pd
from workflow_checks import validate_table, variant_keys, record_qc
import yaml

# Core columns needed by the current Stage-3 filtering checkpoint.
KEY_COLS = ["Chr_Pos", "Variant_REF_ALT"]
TX_SEP = " | "
MISSING_VALUES = {"", ".", "-", "_", "nan", "none", "<na>"}
VEP_URL = "https://rest.ensembl.org/vep/homo_sapiens/region"
VEP_BATCH_SIZE = 100

VEP_CACHE_VERSION = 1

STAGE3_VEP_CACHE = Path(
    "results/intermediate/stage3/vep_after_rarity.json.gz"
)

FINAL_COLUMNS = [
    # Variant identity
    "Chr_Pos",
    "Variant_REF_ALT",

    # Gene / transcript annotation
    "Gene_symbol",
    "Gene_ID",
    "Transcript_ID",
    "MANE_Transcript_ID",
    "Transcript_biotype",
    "Consequences",
    "HGVSc",
    "HGVSp",
    "Gene_context",
    "Regulatory_context",

    # Proband
    "GT",
    "GQ",
    "DP",
    "AD_alt",
    "VAF",

    # Nuclear gnomAD
    "gnomAD_AF",
    "gnomAD_grpmax_AF",
    "gnomAD_het",
    "gnomAD_hom",
    "gnomAD_hemi",

    # Mitochondrial gnomAD
    "gnomAD_mt_AF_hom",
    "gnomAD_mt_AF_het",

    # ClinVar
    "ClinVar_significance",

    # Case / secondary context
    "Case_gene_match",
    "Secondary_gene_match",
    "Secondary_MOI",
    "Secondary_reportable",

    # Distance annotation
    "Splice_context",
    "Nearest_splice_distance",
    "Nearest_splice_transcript",
    "TSS_context",
    "Nearest_TSS_distance",
    "Nearest_TSS_transcript",
    "Nearest_TSS_gene",
    "Regulatory_relevant",

    # Compound het
    "Compound_het_gene",
    "Compound_het_partner",

    # Keep last
    "inheritance_category",
]

REQUIRED_COLUMNS = [
    "Chr_Pos", "Variant_REF_ALT", "FILTER", "GT", "GQ", "DP", "AD_alt", "VAF",
    "gnomAD_AF", "gnomAD_grpmax_AF", "gnomAD_het", "gnomAD_hom", "gnomAD_hemi",
    "gnomAD_mt_AF_hom", "gnomAD_mt_AF_het", "ClinVar_significance",
    "Gene_symbol", "Gene_ID", "Transcript_ID", "Transcript_biotype", "Consequences",
    "Regulatory_context", "MANE_Transcript_ID",
]

# Functional rules
CODING_LIKE_BIOTYPES = {
    "protein_coding", "nonsense_mediated_decay", "non_stop_decay",
    "protein_coding_CDS_not_defined", "protein_coding_LoF",
}
MEANINGFUL_CODING_CONSEQUENCES = {
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
    "splice_region_variant", "splice_donor_region_variant",
    "splice_donor_5th_base_variant", "splice_polypyrimidine_tract_variant",
    "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
    "missense_variant", "inframe_insertion", "inframe_deletion",
    "protein_altering_variant", "coding_sequence_variant","synonymous_variant",
}
SPECIAL_NONCODING_CONSEQUENCES = {("miRNA", "mature_miRNA_variant")}
DIRECT_ENCODE_CLASSES = {"PLS", "pELS"}
ENCODE_ENHANCER_SUPPORT_CLASSES = {"PLS", "pELS"}
PSEUDOGENE_TOKEN = "pseudogene"

# ACMG SF v3.3: heterozygous females are individually reportable for these
# X-linked conditions. Other XL genes require hemizygous/homozygous or a
# biallelic P/LP combination.
XL_SINGLE_HET_REPORTABLE_GENES = {"GLA", "OTC"}
ACMG_MOI_VALUES = {"AD", "AR", "SD", "XL"}

# QC / inheritance thresholds
QC_GQ_MIN = 20
QC_DP_MIN = 10
MOSAIC_DP_MIN = 20
MOSAIC_ALT_DEPTH_MIN = 5
MOSAIC_VAF_MIN = 0.10
MOSAIC_VAF_MAX = 0.35
MOSAIC_AF_MAX = 0.001
MITO_AF_MAX = 0.005
MITO_HETEROPLASMY_MIN = 0.10
GERMLINE_HET_VAF_MIN = 0.35
GERMLINE_HET_VAF_MAX = 0.65

RARITY_THRESHOLDS = {
    # Reserved for the later final heterozygous stage.
    "AD_singleton": {"af_max": 1e-4, "het_lt": 30, "hom_max": 0},
    "AD_case_singleton": {"af_max": 1e-3, "het_lt": 50, "hom_max": 0},
    # Used now to build the preliminary heterozygous pool.
    "AR_compound_generic": {"af_max": 1e-3, "het_lt": 500, "hom_max": 2},
    "AR_compound_targeted": {"af_max": 1e-2, "het_lt": None, "hom_max": 2},
    "AR_homozygous": {"af_max": 1e-3, "het_lt": 500, "hom_max": 2},
}

AUTOSOMES = {str(i) for i in range(1, 23)}
PAR1 = (10001, 2781479)  # GRCh38, 1-based inclusive
PAR2 = (155701383, 156030895)
SPLICE_DISTANCE_MAX = 50
TSS_DISTANCE_MAX = 10_000

GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')
TX_ID_RE = re.compile(r'transcript_id "([^"]+)"')
ENCODE_CLASS_RE = re.compile(
    r"ENCODE:\s*(?:[^:;|]+\s*:\s*)?(PLS|pELS|dELS)\b",
    re.IGNORECASE,
)


def is_missing(value):
    if pd.isna(value):
        return True
    return str(value).strip().lower() in MISSING_VALUES


def split_tx_field(value):
    """Split aligned transcript fields while preserving missing placeholders."""
    if pd.isna(value):
        return []
    return [part.strip() for part in str(value).split(TX_SEP)]


def consequence_terms(value):
    if is_missing(value):
        return set()
    return {term.strip() for term in str(value).split("&") if term.strip()}


def strip_version(value):
    if is_missing(value):
        return ""
    return str(value).strip().split(".")[0]


def variant_genes(value):
    """Return unique annotated gene symbols for one collapsed variant."""
    return {
        gene
        for gene in split_tx_field(value)
        if not is_missing(gene)
    }

def parse_gt(gt):
    """Return (het, hom_alt, haploid_alt, has_alt); partial calls stay unknown."""
    if is_missing(gt):
        return False, False, False, False

    alleles = re.split(r"[/|]", str(gt).strip())
    if not alleles or any(a == "." for a in alleles):
        return False, False, False, False

    try:
        alleles = [int(a) for a in alleles]
    except ValueError:
        return False, False, False, False

    has_alt = any(a > 0 for a in alleles)
    if len(alleles) == 1:
        return False, False, has_alt, has_alt

    hom_alt = all(a > 0 for a in alleles) and len(set(alleles)) == 1
    het = has_alt and len(set(alleles)) > 1
    return het, hom_alt, False, has_alt


def add_genotype_flags(df):
    gt = df["GT"].map(parse_gt)
    df[["_het", "_hom_alt", "_haploid_alt", "_has_alt"]] = pd.DataFrame(
        gt.tolist(),
        index=df.index,
        columns=["_het", "_hom_alt", "_haploid_alt", "_has_alt"],
    )
    return df


def clinvar_plp_mask(series):
    """Accept only a pure P, LP or combined P/LP main classification."""
    return series.map(is_main_plp)


def load_collapsed(path):
    """Load Stage-3 input and physically remove QC failures before biology."""
    df = pd.read_parquet(path, columns=REQUIRED_COLUMNS)

    validate_table(df, "Stage 3 collapsed input", REQUIRED_COLUMNS, aligned=True)
    if df[KEY_COLS].duplicated().any():
        raise ValueError("Collapsed parquet is not one-row-per-variant.")

    # The collapsed table already contains the reference patch's keyed MANE
    # values. Carry them directly rather than reselecting transcripts.
    df["MANE_Transcript_ID"] = df["MANE_Transcript_ID"].fillna("").astype(str).str.strip()

    # QC is the first filter. P/LP or any biological annotation cannot rescue it.
    df[["GQ", "DP"]] = df[["GQ", "DP"]].apply(pd.to_numeric, errors="coerce")
    qc = (
        df["FILTER"].fillna("").astype(str).str.upper().eq("PASS")
        & df["GQ"].ge(QC_GQ_MIN)
        & df["DP"].ge(QC_DP_MIN)
    )
    before = len(df)
    df = df.loc[qc].copy()

    df["_key"] = df["Chr_Pos"].astype(str) + "|" + df["Variant_REF_ALT"].astype(str)
    chrom_pos = df["Chr_Pos"].astype(str).str.rsplit(":", n=1, expand=True)
    df["_chrom"] = chrom_pos[0].str.removeprefix("chr").replace({"M": "MT"})
    df["_pos"] = pd.to_numeric(chrom_pos[1], errors="coerce")

    numeric = [
        "AD_alt", "VAF", "gnomAD_AF", "gnomAD_grpmax_AF", "gnomAD_het",
        "gnomAD_hom", "gnomAD_hemi", "gnomAD_mt_AF_hom", "gnomAD_mt_AF_het",
    ]
    df[numeric] = df[numeric].apply(pd.to_numeric, errors="coerce")

    df["_filter_af"] = df[["gnomAD_AF", "gnomAD_grpmax_AF"]].max(axis=1, skipna=True)
    df["_mt_filter_af"] = df[["gnomAD_mt_AF_hom", "gnomAD_mt_AF_het"]].max(axis=1, skipna=True)
    df["_plp"] = clinvar_plp_mask(df["ClinVar_significance"])
    df["_is_par"] = df["_chrom"].eq("X") & (
        df["_pos"].between(*PAR1, inclusive="both")
        | df["_pos"].between(*PAR2, inclusive="both")
    )

    print("\n=== Stage 3 input / QC ===")
    print(f"Collapsed variants: {before:,}")
    print(f"QC passed:          {len(df):,}")
    print(f"QC removed:         {before - len(df):,}")
    return df, before

def classify_function(
    gene_value,
    tx_value,
    biotype_value,
    consequence_value,
    regulatory_value,
):
    """
    Return:
        direct_keep,
        needs_distance,
        intronic,
        dELS,
        has_regulatory

    Rules:
      - Coding / strong splice / mature miRNA -> DIRECT_KEEP
      - ENCODE PLS / pELS -> DIRECT_KEEP
      - Ensembl promoter -> DIRECT_KEEP
      - dELS -> NEEDS_DISTANCE, never DIRECT_KEEP merely because an
        Ensembl enhancer annotation is also present
      - intronic -> NEEDS_DISTANCE unless another direct rule already applies
    """

    genes = split_tx_field(gene_value)
    transcripts = split_tx_field(tx_value)
    biotypes = split_tx_field(biotype_value)
    consequences = split_tx_field(consequence_value)

    if len(biotypes) != len(consequences):
        raise ValueError(
            "Transcript_biotype and Consequences lost alignment."
        )

    meaningful = False
    intronic = False

    for biotype, consequence in zip(biotypes, consequences):
        terms = consequence_terms(consequence)

        intronic |= "intron_variant" in terms

        if (
            biotype in CODING_LIKE_BIOTYPES
            and terms & MEANINGFUL_CODING_CONSEQUENCES
        ):
            meaningful = True

        elif any(
            biotype == special_biotype
            and special_consequence in terms
            for special_biotype, special_consequence
            in SPECIAL_NONCODING_CONSEQUENCES
        ):
            meaningful = True

    has_gene = any(
        not is_missing(x)
        for x in genes
    )

    has_transcript = any(
        not is_missing(x)
        for x in transcripts
    )

    intergenic = not has_gene and not has_transcript

    nonmissing_biotypes = [
        x
        for x in biotypes
        if not is_missing(x)
    ]

    pseudogene_only = (
        bool(nonmissing_biotypes)
        and all(
            PSEUDOGENE_TOKEN in x.lower()
            for x in nonmissing_biotypes
        )
    )

    has_regulatory = not is_missing(regulatory_value)

    encode_classes = set()
    has_dels = False
    strong_regulatory = False

    if has_regulatory:

        text = str(regulatory_value)
        low = text.lower()

        encode_classes = {
            {
                "pls": "PLS",
                "pels": "pELS",
                "dels": "dELS",
            }[match.lower()]
            for match in ENCODE_CLASS_RE.findall(text)
        }

        has_dels = "dELS" in encode_classes

        ensembl_promoter = (
            "ensembl:" in low
            and "promoter" in low
        )

        # IMPORTANT:
        # dELS is NOT direct keep.
        #
        # If an annotation contains:
        #   ENCODE dELS + Ensembl enhancer
        # it still requires the TSS distance rule.
        #
        # PLS/pELS remain direct because the ENCODE class itself
        # is considered sufficiently strong.
        strong_regulatory = (
            bool(encode_classes & DIRECT_ENCODE_CLASSES)
            or ensembl_promoter
        )

    direct_keep = (
        meaningful
        or strong_regulatory
    )

    needs_distance = (
        not direct_keep
        and (intronic or has_dels)
        and not intergenic
        and not pseudogene_only
    )

    return (
        direct_keep,
        needs_distance,
        intronic,
        has_dels,
        has_regulatory,
    )

def run_functional_classification(df):
    """
    Split QC-passing variants into functional keep categories.

    A variant survives this stage if it is DIRECT_KEEP, NEEDS_DISTANCE, or
    ClinVar P/LP. P/LP is a separate functional keep category, not a QC rescue.
    Everything else is DROP and is removed here.
    """
    flag_cols = ["_direct_keep", "_needs_distance", "_intronic", "_dels", "_has_regulatory"]
    values = [
        classify_function(gene, tx, biotype, consequence, regulatory)
        for gene, tx, biotype, consequence, regulatory in zip(
            df["Gene_symbol"], df["Transcript_ID"], df["Transcript_biotype"],
            df["Consequences"], df["Regulatory_context"],
        )
    ]
    df[flag_cols] = pd.DataFrame(values, index=df.index, columns=flag_cols)
    # P/LP is an independent functional keep category after QC. Secondary
    # classification is added later and never removes the underlying P/LP fact.

    functional_keep = df["_direct_keep"] | df["_needs_distance"] | df["_plp"]

    print("\n=== Functional classification ===")
    print(f"QC-passing input:   {len(df):,}")
    print(f"DIRECT_KEEP:        {int(df['_direct_keep'].sum()):,}")
    print(f"NEEDS_DISTANCE:     {int(df['_needs_distance'].sum()):,}")
    print(f"P/LP keep category: {int(df['_plp'].sum()):,}")
    print(f"DROP removed:       {int((~functional_keep).sum()):,}")
    return df.loc[functional_keep].copy()

def load_case_context(path):
    with open(path) as f:
        context = yaml.safe_load(f) or {}

    sex = str(context.get("sex", "unknown")).strip().lower()

    if sex not in {"male", "female"}:
        raise ValueError(
            f"case_context sex must be 'male' or 'female', got: {sex!r}"
        )

    case_genes = {
        str(g).strip()
        for g in context.get("case_genes", [])
        if str(g).strip()
    }

    # Optional case-specific override.
    # Applies to recessive homozygote and X-linked hemizygote limits, not AD.
    hom_max_override = (
        context
        .get("hom_max_override")
    )

    if hom_max_override is not None:
        if isinstance(hom_max_override, bool) or not isinstance(hom_max_override, (int, str)):
            raise ValueError("hom_max_override must be a nonnegative integer")
        try:
            hom_max_override = int(hom_max_override)
        except ValueError:
            raise ValueError("hom_max_override must be a nonnegative integer") from None
        if hom_max_override < 0:
            raise ValueError(
                "hom_max_override must be >= 0"
            )

    return sex, case_genes, hom_max_override

def add_case_and_genotype_context(df, sex, case_genes):
    """Add genotype plus case/P-LP-gene flags to the surviving functional pool."""

    plp_genes = {
        gene
        for value in df.loc[df["_plp"], "Gene_symbol"]
        for gene in split_tx_field(value)
        if not is_missing(gene)
    }

    def gene_flags(value):
        genes = {g for g in split_tx_field(value) if not is_missing(g)}
        return bool(genes & case_genes), bool(genes & plp_genes)

    flags = df["Gene_symbol"].map(gene_flags)
    df[["_case_gene", "_plp_gene"]] = pd.DataFrame(
        flags.tolist(), index=df.index, columns=["_case_gene", "_plp_gene"]
    )
    df["_second_hit"] = df["_plp_gene"] & ~df["_plp"]

    print("\n=== Case / P-LP gene context ===")
    print(f"Proband sex:     {sex}")
    print(f"Case genes:      {len(case_genes):,}")
    print(f"Genes with P/LP: {len(plp_genes):,}")
    return df


def load_acmg_sf(path):
    """Return gene -> set of ACMG SF inheritance modes."""
    acmg = pd.read_csv(path, dtype=str)
    required = {"Gene Symbol", "MOI", "Reportable as SF"}

    if not required.issubset(acmg.columns):
        acmg = pd.read_csv(path, dtype=str, skiprows=1)
    if not required.issubset(acmg.columns):
        raise ValueError(f"ACMG file must contain {sorted(required)}")

    acmg = acmg[
        acmg["Reportable as SF"]
        .fillna("")
        .str.strip()
        .isin(["Yes", "See Guidance"])
    ]

    gene_moi = defaultdict(set)
    unknown_modes = set()

    for gene, moi in zip(acmg["Gene Symbol"], acmg["MOI"]):
        if is_missing(gene) or is_missing(moi):
            continue

        modes = {
            token
            for token in re.split(r"[/,;|\s]+", str(moi).upper())
            if token
        }
        unknown_modes.update(modes - ACMG_MOI_VALUES)
        gene_moi[str(gene).strip()].update(modes & ACMG_MOI_VALUES)

    if unknown_modes:
        raise ValueError(
            "Unexpected ACMG MOI value(s): " + ", ".join(sorted(unknown_modes))
        )

    return gene_moi


def add_secondary_context(df, sex, gene_moi):
    """
    Annotate ACMG Secondary-gene membership and automated reportability.

    Secondary_gene_match is annotation only. Secondary_reportable requires P/LP
    plus an ACMG-compatible allelic state. For recessive disease, two P/LP
    heterozygous variants at different loci qualify as a biallelic candidate;
    phase is not assumed. For X-linked disease, a single heterozygous female is
    individually reportable only for GLA and OTC under ACMG SF v3.3.
    """
    # Precompute heterozygous P/LP alleles per gene so AR and XL biallelic
    # Secondary logic can be evaluated without repeated whole-dataframe scans.
    plp_het_by_gene = defaultdict(dict)
    for key, locus, gene_value, is_plp, is_het in zip(
        df["_key"],
        df["Chr_Pos"],
        df["Gene_symbol"],
        df["_plp"],
        df["_het"],
    ):
        if not (is_plp and is_het):
            continue
        for gene in variant_genes(gene_value):
            plp_het_by_gene[gene][key] = str(locus)

    def has_distinct_plp_partner(gene, key, locus):
        return any(
            partner_key != key and partner_locus != str(locus)
            for partner_key, partner_locus in plp_het_by_gene.get(gene, {}).items()
        )

    secondary_gene = []
    secondary_moi = []
    secondary_reportable = []

    for _, row in df.iterrows():
        genes = variant_genes(row["Gene_symbol"])
        matched_mois = set()
        reportable = False

        for gene in genes:
            mois = gene_moi.get(gene, set())
            if not mois:
                continue

            matched_mois.update(mois)
            if not row["_plp"]:
                continue

            # Autosomal dominant / semidominant: one P/LP ALT allele is enough.
            if mois & {"AD", "SD"} and (row["_het"] or row["_hom_alt"]):
                reportable = True

            # Autosomal recessive: homozygous P/LP or two P/LP heterozygous
            # variants at distinct genomic loci. Phase is not claimed here.
            if "AR" in mois:
                if row["_hom_alt"]:
                    reportable = True
                elif row["_het"] and has_distinct_plp_partner(
                    gene, row["_key"], row["Chr_Pos"]
                ):
                    reportable = True

            # X-linked ACMG SF v3.3.
            if "XL" in mois and row["_chrom"] == "X" and not row["_is_par"]:
                if sex == "male":
                    if row["_hom_alt"] or row["_haploid_alt"]:
                        reportable = True
                else:
                    if row["_hom_alt"]:
                        reportable = True
                    elif row["_het"] and gene in XL_SINGLE_HET_REPORTABLE_GENES:
                        reportable = True
                    elif row["_het"] and has_distinct_plp_partner(
                        gene, row["_key"], row["Chr_Pos"]
                    ):
                        reportable = True

        secondary_gene.append(bool(matched_mois))
        secondary_moi.append(TX_SEP.join(sorted(matched_mois)))
        secondary_reportable.append(reportable)

    df["_secondary_gene"] = secondary_gene
    df["_secondary_moi"] = secondary_moi
    df["_secondary_reportable"] = secondary_reportable

    print("\n=== ACMG Secondary ===")
    print(f"Variants in ACMG genes: {int(df['_secondary_gene'].sum()):,}")
    print(f"Reportable P/LP:        {int(df['_secondary_reportable'].sum()):,}")
    return df


def run_initial_rarity_filters(df, sex, hom_max_override=None):
    """
    Build the preliminary genotype/rarity pool.

    P/LP variants already passed QC and were retained as their own functional
    keep category, so population rarity does not remove them here. Rows that fail
    every route are removed before the more expensive distance annotation.
    """

    hom_max = (
        hom_max_override
        if hom_max_override is not None
        else RARITY_THRESHOLDS["AR_homozygous"]["hom_max"]
    )
    x_nonpar = df["_chrom"].eq("X") & ~df["_is_par"]
    diploid_region = (
        df["_chrom"].isin(AUTOSOMES) | df["_chrom"].eq("X")
        if sex == "female"
        else df["_chrom"].isin(AUTOSOMES) | df["_is_par"]
    )
    # Keep this explicit for later compound-het pairing. A heterozygous P/LP
    # variant may bypass rarity, but it must still be in a biologically diploid
    # region to participate in a compound-heterozygous pair.
    df["_diploid_region"] = diploid_region

    hom_rule = RARITY_THRESHOLDS["AR_homozygous"]
    hom_af_ok = df["_filter_af"].le(hom_rule["af_max"]) | df["_filter_af"].isna()
    het_count_ok = df["gnomAD_het"].lt(hom_rule["het_lt"]) | df["gnomAD_het"].isna()
    hom_count_ok = df["gnomAD_hom"].le(hom_max) | df["gnomAD_hom"].isna()

    df["_hom"] = df["_hom_alt"] & diploid_region & hom_af_ok & het_count_ok & hom_count_ok
    df["_x_hemi"] = False
    if sex == "male":
        hemi_ok = df["gnomAD_hemi"].le(hom_max) | df["gnomAD_hemi"].isna()
        df["_x_hemi"] = (
                x_nonpar & (df["_hom_alt"] | df["_haploid_alt"])
                & hom_af_ok & het_count_ok & hemi_ok
        )

    df["_mosaic"] = (
            ~df["_chrom"].eq("MT") & df["_has_alt"]
            & df["DP"].ge(MOSAIC_DP_MIN) & df["AD_alt"].ge(MOSAIC_ALT_DEPTH_MIN)
            & df["VAF"].between(MOSAIC_VAF_MIN, MOSAIC_VAF_MAX, inclusive="both")
            & (df["_filter_af"].le(MOSAIC_AF_MAX) | df["_filter_af"].isna())
    )
    df["_mito"] = (
            df["_chrom"].eq("MT") & df["_has_alt"]
            & df["VAF"].ge(MITO_HETEROPLASMY_MIN)
            & (df["_mt_filter_af"].le(MITO_AF_MAX) | df["_mt_filter_af"].isna())
    )

    base_het = (
        df["_het"]
        & diploid_region
        & ~df["_plp"]
        & ~df["_mosaic"]
        & df["VAF"].between(
            GERMLINE_HET_VAF_MIN,
            GERMLINE_HET_VAF_MAX,
            inclusive="both",
        )
    )
    generic_rule = RARITY_THRESHOLDS["AR_compound_generic"]
    targeted_rule = RARITY_THRESHOLDS["AR_compound_targeted"]

    generic = (
            base_het & ~df["_case_gene"] & ~df["_plp_gene"]
            & (df["_filter_af"].le(generic_rule["af_max"]) | df["_filter_af"].isna())
            & (df["gnomAD_het"].lt(generic_rule["het_lt"]) | df["gnomAD_het"].isna())
            & (df["gnomAD_hom"].le(hom_max) | df["gnomAD_hom"].isna())
    )
    targeted = (
            base_het & (df["_case_gene"] | df["_plp_gene"])
            & (df["_filter_af"].le(targeted_rule["af_max"]) | df["_filter_af"].isna())
            & (df["gnomAD_hom"].le(hom_max) | df["gnomAD_hom"].isna())
    )
    df["_het_eligible"] = generic | targeted

    keep = (
            df["_plp"]
            | df["_hom"]
            | df["_x_hemi"]
            | df["_mosaic"]
            | df["_mito"]
            | df["_het_eligible"]
    )

    print("\n=== Initial genotype + rarity ===")
    print(f"Homozygous:       {int(df['_hom'].sum()):,}")
    print(f"X hemizygous:     {int(df['_x_hemi'].sum()):,}")
    print(f"Mosaic:           {int(df['_mosaic'].sum()):,}")
    print(f"Mitochondrial:    {int(df['_mito'].sum()):,}")
    print(f"Het generic:      {int(generic.sum()):,}")
    print(f"Het targeted:     {int(targeted.sum()):,}")
    print(f"P/LP retained:    {int(df['_plp'].sum()):,}")
    print(f"Retained after rarity: {int(keep.sum()):,}")
    return df.loc[keep].copy()

def call_vep(variants):

    params = {
        "hgvs": 1,
        "numbers": 1,
        "transcript_version": 1,
        "distance": TSS_DISTANCE_MAX,
        "domains": 1,

        # Already request Stage 4 evidence now.
        # We do not need to use all of it in Stage 3.
        "CADD": 1,
        "REVEL": 1,
        "AlphaMissense": 1,
        "LoF": 1,
        "NMD": 1,
    }

    for attempt in range(1, 6):

        response = requests.post(
            VEP_URL,
            params=params,
            json={"variants": variants},
            headers={
                "Content-Type": "application/json",
                "Accept": "application/json",
            },
            timeout=300,
        )

        if response.ok:
            return response.json()

        print(
            f"VEP attempt {attempt}/5 failed: "
            f"HTTP {response.status_code}"
        )

        if response.status_code == 429:
            wait = int(
                response.headers.get(
                    "Retry-After",
                    60,
                )
            )
            print(
                f"Rate limited. Waiting {wait}s..."
            )
            time.sleep(wait)
            continue

        if response.status_code >= 500:
            wait = 15 * attempt
            print(
                f"Server error. Waiting {wait}s..."
            )
            time.sleep(wait)
            continue

        response.raise_for_status()

    raise RuntimeError(
        "VEP failed after 5 attempts."
    )

def save_vep_cache(cache_path, vep_input, results, complete=False):
    cache_path = Path(cache_path)
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    payload = {
        "version": VEP_CACHE_VERSION,
        "vep_input": vep_input,
        "complete": complete,
        "results": results,
    }

    tmp_path = cache_path.with_suffix(
        cache_path.suffix + ".tmp"
    )

    with gzip.open(tmp_path, "wt") as f:
        json.dump(payload, f)

    tmp_path.replace(cache_path)

def vep_input_key(value):
    """
    Stable VEP variant key independent of the temporary S3_* ID.
    """

    parts = str(value).split()

    if len(parts) < 3:
        return None

    return tuple(
        parts[:2] + parts[3:]
    )

def load_or_run_stage3_vep(df, cache_path):

    print("\n" + "=" * 60)
    print("STAGE 3 — VEP AFTER RARITY")
    print("=" * 60)

    vep_input = []
    id_to_variant = {}

    for n, (_, row) in enumerate(
        df.iterrows()
    ):

        chrom = str(row["_chrom"])
        pos = int(row["_pos"])

        # mtDNA will be handled separately later.
        if chrom == "MT":
            continue

        ref, alt = re.split(
            r"[>_/]",
            str(row["Variant_REF_ALT"]),
            maxsplit=1,
        )

        variant_id = f"S3_{n}"

        allowed_txs = {
            strip_version(tx)
            for tx in split_tx_field(
                row["Transcript_ID"]
            )
            if not is_missing(tx)
        }

        id_to_variant[variant_id] = {
            "Chr_Pos": row["Chr_Pos"],
            "Variant_REF_ALT":
                row["Variant_REF_ALT"],
            "chrom": chrom,
            "pos": pos,
            "transcripts": allowed_txs,
        }

        vep_input.append(
            f"{chrom} {pos} {variant_id} "
            f"{ref} {alt} . . ."
        )

    print(
        f"Variants sent to nuclear VEP: "
        f"{len(vep_input):,}"
    )

    # --------------------------------------------------------
    # Try cache
    # --------------------------------------------------------

    cache_path = Path(cache_path)

    # ========================================================
    # Load existing complete or partial cache
    # ========================================================

    vep_results = []

    if cache_path.exists():

        print(
            f"Found existing VEP cache: "
            f"{cache_path}"
        )

        try:
            with gzip.open(
                cache_path,
                "rt"
            ) as f:
                cache_data = json.load(f)

            cached_results = (
                cache_data.get("results", [])
            )

            if (
                cache_data.get("version")
                != VEP_CACHE_VERSION
                or not isinstance(
                    cached_results,
                    list,
                )
            ):
                print(
                    "Existing VEP cache version/results "
                    "are incompatible."
                )
                print(
                    "Starting a new VEP run."
                )

            elif (
                cache_data.get("vep_input")
                == vep_input
            ):
                # Exact same input
                vep_results = cached_results

                print(
                    f"Loaded VEP cache: "
                    f"{len(vep_results):,} results"
                )

            else:
                # Input changed slightly:
                # reuse overlapping variants
                print(
                    "VEP input changed — checking "
                    "for reusable cached variants."
                )

                cached_by_key = {}

                for result in cached_results:

                    key = vep_input_key(
                        result.get(
                            "input",
                            "",
                        )
                    )

                    if key is not None:
                        cached_by_key[key] = result


                for current_input in vep_input:

                    key = vep_input_key(
                        current_input
                    )

                    cached_result = (
                        cached_by_key.get(key)
                    )

                    if cached_result is None:
                        continue

                    reused = dict(
                        cached_result
                    )

                    # Replace old S3_* ID with
                    # the ID from the current run
                    reused["input"] = (
                        current_input
                    )

                    vep_results.append(
                        reused
                    )

                print(
                    f"Reused overlapping VEP results: "
                    f"{len(vep_results):,}"
                )

        except Exception as e:

            print(
                f"Could not read existing VEP cache: "
                f"{e}"
            )

            print(
                "Starting a new VEP run."
            )

    # ========================================================
    # Determine which variants are already completed
    # ========================================================

    expected_ids = {
        variant.split()[2]
        for variant in vep_input
    }

    completed_ids = set()

    for result in vep_results:

        parts = str(
            result.get("input", "")
        ).split()

        if len(parts) >= 3:
            completed_ids.add(parts[2])


    missing_input = [
        variant
        for variant in vep_input
        if variant.split()[2]
        not in completed_ids
    ]


    print("\nVEP progress")
    print("-" * 60)
    print(
        f"Expected variants:  "
        f"{len(expected_ids):,}"
    )
    print(
        f"Already completed: "
        f"{len(completed_ids):,}"
    )
    print(
        f"Remaining:         "
        f"{len(missing_input):,}"
    )


    # ========================================================
    # Run only missing batches
    # ========================================================

    if missing_input:

        total_batches = (
            len(missing_input)
            + VEP_BATCH_SIZE
            - 1
        ) // VEP_BATCH_SIZE

        for start in range(
            0,
            len(missing_input),
            VEP_BATCH_SIZE,
        ):

            batch = missing_input[
                start:
                start + VEP_BATCH_SIZE
            ]

            batch_number = (
                start // VEP_BATCH_SIZE
            ) + 1

            print(
                f"\nVEP batch "
                f"{batch_number}/{total_batches} "
                f"({len(batch)} variants)"
            )

            batch_results = call_vep(batch)

            vep_results.extend(
                batch_results
            )

            # Save checkpoint after every successful batch
            save_vep_cache(
                cache_path,
                vep_input,
                vep_results,
                complete=False,
            )

            print(
                f"Checkpoint saved: "
                f"{len(vep_results):,} "
                f"VEP results"
            )

            if start + VEP_BATCH_SIZE < len(
                missing_input
            ):
                time.sleep(1)

    else:
        print(
            "VEP cache already contains "
            "all expected variants."
        )


    # ========================================================
    # Final validation
    # ========================================================

    returned_ids = set()

    for result in vep_results:

        parts = str(
            result.get("input", "")
        ).split()

        if len(parts) >= 3:
            returned_ids.add(parts[2])


    missing_ids = (
        expected_ids
        - returned_ids
    )

    extra_ids = (
        returned_ids
        - expected_ids
    )

    if missing_ids or extra_ids:

        raise RuntimeError(
            "VEP ID validation failed: "
            f"missing={len(missing_ids):,}, "
            f"extra={len(extra_ids):,}"
        )


    # ========================================================
    # Mark cache as complete
    # ========================================================

    save_vep_cache(
        cache_path,
        vep_input,
        vep_results,
        complete=True,
    )

    print(
        f"\nVEP complete: "
        f"{len(returned_ids):,} / "
        f"{len(expected_ids):,} variants"
    )

    return vep_results, id_to_variant

def load_distance_resources(
    gtf_path,
    transcript_ids,
):
    """
    Load transcript-specific TSSs and splice boundaries
    only for transcripts present in current candidates.
    """

    transcript_ids = {
        strip_version(tx)
        for tx in transcript_ids
        if tx
    }

    tss_by_tx = {}
    exons_by_tx = defaultdict(list)

    opener = (
        gzip.open
        if str(gtf_path).endswith(".gz")
        else open
    )

    with opener(gtf_path, "rt") as handle:

        for line in handle:

            if line.startswith("#"):
                continue

            parts = line.rstrip().split("\t")

            if len(parts) != 9:
                continue

            (
                chrom,
                _,
                feature,
                start,
                end,
                _,
                strand,
                _,
                attrs,
            ) = parts

            if feature not in {
                "transcript",
                "exon",
            }:
                continue

            match = TX_ID_RE.search(attrs)

            if not match:
                continue

            tx = strip_version(
                match.group(1)
            )

            if tx not in transcript_ids:
                continue

            chrom = chrom.removeprefix(
                "chr"
            )

            if chrom == "M":
                chrom = "MT"

            if feature == "transcript":

                tss = (
                    int(start)
                    if strand == "+"
                    else int(end)
                )

                tss_by_tx[tx] = (
                    chrom,
                    tss,
                    strand,
                )

            else:

                exons_by_tx[tx].append(
                    (
                        chrom,
                        int(start),
                        int(end),
                    )
                )

    splice_by_tx = {}

    for tx, exons in exons_by_tx.items():

        exons = sorted(
            exons,
            key=lambda x: (
                x[0],
                x[1],
                x[2],
            ),
        )

        if len(exons) < 2:
            continue

        sites = []

        for i, (
            chrom,
            start,
            end,
        ) in enumerate(exons):

            if i > 0:
                sites.append(
                    (chrom, start)
                )

            if i < len(exons) - 1:
                sites.append(
                    (chrom, end)
                )

        splice_by_tx[tx] = sites

    return tss_by_tx, splice_by_tx

def build_transcript_context(
    df,
    vep_results,
    id_to_variant,
    gtf_path,
):
    """
    Build temporary variant × transcript context from Stage 3 VEP results
    and add transcript-specific splice/TSS distances.

    Returns one row per variant × transcript.
    """

    print("\n=== Building transcript context ===")

    # --------------------------------------------------------
    # Build Stage 3 variant × transcript map
    # --------------------------------------------------------

    expected = {}
    allowed_transcripts = set()

    for _, row in df.iterrows():

        chr_pos = row["Chr_Pos"]
        ref_alt = row["Variant_REF_ALT"]

        txs = split_tx_field(
            row.get("Transcript_ID")
        )

        consequences = split_tx_field(
            row.get("Consequences")
        )

        if len(txs) != len(consequences):
            raise RuntimeError(
                "Transcript/consequence alignment failed for "
                f"{chr_pos} {ref_alt}: "
                f"{len(txs)} transcripts vs "
                f"{len(consequences)} consequences"
            )

        for tx_raw, consequence in zip(
            txs,
            consequences,
        ):

            tx = strip_version(tx_raw)

            if not tx:
                continue

            key = (
                chr_pos,
                ref_alt,
                tx,
            )

            expected[key] = {
                "Transcript_ID": tx_raw,
                "Consequence": consequence,
            }

            allowed_transcripts.add(tx)


    print(
        f"Stage 3 transcript associations: "
        f"{len(expected):,}"
    )

    print(
        f"Unique transcripts: "
        f"{len(allowed_transcripts):,}"
    )


    # --------------------------------------------------------
    # Parse VEP transcript results
    # --------------------------------------------------------

    rows = []
    observed = set()

    for result in vep_results:

        input_parts = str(
            result.get("input", "")
        ).split()

        if len(input_parts) < 3:
            continue

        variant_id = input_parts[2]

        if variant_id not in id_to_variant:
            continue

        variant_info = id_to_variant[
            variant_id
        ]

        chr_pos = variant_info["Chr_Pos"]
        ref_alt = variant_info["Variant_REF_ALT"]
        chrom = str(
            variant_info["chrom"]
        )
        pos = int(
            variant_info["pos"]
        )

        for tc in result.get(
            "transcript_consequences",
            [],
        ):

            tx = strip_version(
                tc.get("transcript_id")
            )

            if not tx:
                continue

            key = (
                chr_pos,
                ref_alt,
                tx,
            )

            # Keep only transcript associations
            # already present in Stage 3
            if key not in expected:
                continue

            observed.add(key)

            # --------------------------------------------
            # Gene / transcript context
            # --------------------------------------------

            gene_context = pd.NA

            if tc.get("exon"):
                gene_context = (
                    f"exon {tc['exon']}"
                )

            elif tc.get("intron"):
                gene_context = (
                    f"intron {tc['intron']}"
                )

            else:
                terms = set(
                    tc.get(
                        "consequence_terms",
                        [],
                    )
                )

                distance = tc.get(
                    "distance"
                )

                if "5_prime_UTR_variant" in terms:
                    gene_context = "5_prime_UTR"

                elif "3_prime_UTR_variant" in terms:
                    gene_context = "3_prime_UTR"

                elif "upstream_gene_variant" in terms:
                    if distance is not None:
                        gene_context = (
                            f"upstream {distance} bp"
                        )
                    else:
                        gene_context = "upstream"

                elif "downstream_gene_variant" in terms:
                    if distance is not None:
                        gene_context = (
                            f"downstream {distance} bp"
                        )
                    else:
                        gene_context = "downstream"


            rows.append({
                "Chr_Pos":
                    chr_pos,

                "Variant_REF_ALT":
                    ref_alt,

                "Transcript_ID":
                    expected[key][
                        "Transcript_ID"
                    ],

                "_Transcript_base":
                    tx,

                "Gene_ID":
                    tc.get("gene_id"),

                "Consequence":
                    expected[key][
                        "Consequence"
                    ],

                "HGVSc":
                    tc.get("hgvsc"),

                "HGVSp":
                    tc.get("hgvsp"),

                "Gene_context":
                    gene_context,

                "_chrom":
                    chrom,

                "_pos":
                    pos,
            })


    # --------------------------------------------------------
    # Validate VEP transcript recovery
    # --------------------------------------------------------

    missing = (
        set(expected)
        - observed
    )

    extra = (
        observed
        - set(expected)
    )

    print(
        f"Parsed transcript associations: "
        f"{len(observed):,}"
    )

    if missing:
        print(
            f"Missing transcript associations: "
            f"{len(missing):,}"
        )

        for key in list(missing)[:10]:
            print(
                "  missing:",
                key,
            )

    if extra:
        print(
            f"Unexpected transcript associations: "
            f"{len(extra):,}"
        )

    if missing or extra:
        raise RuntimeError(
            "Stage 3 / VEP transcript "
            "association mismatch."
        )


    tx_context = pd.DataFrame(rows)

    if tx_context.empty:
        raise RuntimeError(
            "Transcript context is empty."
        )


    # --------------------------------------------------------
    # Load transcript-specific GENCODE distances
    # --------------------------------------------------------

    tss_by_tx, splice_by_tx = (
        load_distance_resources(
            gtf_path,
            allowed_transcripts,
        )
    )


    # --------------------------------------------------------
    # Normalize chromosomes for comparison
    # --------------------------------------------------------

    def norm_chrom(value):
        value = str(value)

        if value.startswith("chr"):
            value = value[3:]

        if value == "M":
            value = "MT"

        return value


    # --------------------------------------------------------
    # Calculate transcript-specific distances
    # --------------------------------------------------------

    splice_distances = []
    tss_distances = []

    for _, row in tx_context.iterrows():

        tx = row["_Transcript_base"]

        chrom = norm_chrom(
            row["_chrom"]
        )

        pos = int(
            row["_pos"]
        )


        # --------------------------------------------
        # TSS distance
        # --------------------------------------------

        tss_distance = pd.NA

        tss_info = tss_by_tx.get(tx)

        if tss_info is not None:

            (
                tss_chrom,
                tss_pos,
                strand,
            ) = tss_info

            if (
                norm_chrom(tss_chrom)
                == chrom
            ):

                if strand == "+":
                    tss_distance = (
                        pos - int(tss_pos)
                    )
                else:
                    tss_distance = (
                        int(tss_pos) - pos
                    )


        # --------------------------------------------
        # Splice distance
        # --------------------------------------------

        splice_distance = pd.NA

        splice_sites = (
            splice_by_tx.get(
                tx,
                [],
            )
        )

        candidate_distances = []

        for splice_info in splice_sites:

            splice_chrom = splice_info[0]
            splice_pos = int(
                splice_info[1]
            )

            if (
                norm_chrom(splice_chrom)
                != chrom
            ):
                continue

            # Use transcript strand when available.
            # Same sign convention as TSS distance.
            strand = (
                tss_info[2]
                if tss_info is not None
                else "+"
            )

            if strand == "+":
                distance = (
                    pos - splice_pos
                )
            else:
                distance = (
                    splice_pos - pos
                )

            candidate_distances.append(
                distance
            )


        if candidate_distances:

            splice_distance = min(
                candidate_distances,
                key=abs,
            )


        splice_distances.append(
            splice_distance
        )

        tss_distances.append(
            tss_distance
        )


    tx_context[
        "Splice_distance"
    ] = splice_distances

    tx_context[
        "TSS_distance"
    ] = tss_distances


    # --------------------------------------------------------
    # Final cleanup / summary
    # --------------------------------------------------------

    tx_context = tx_context.drop(
        columns=[
            "_Transcript_base",
            "_chrom",
            "_pos",
        ]
    )


    print(
        f"Transcript context rows: "
        f"{len(tx_context):,}"
    )

    print(
        "Splice distances available: "
        f"{tx_context['Splice_distance'].notna().sum():,}"
    )

    print(
        "TSS distances available: "
        f"{tx_context['TSS_distance'].notna().sum():,}"
    )

    print(
        "Transcript context: PASS"
    )

    return tx_context

def add_distance_annotations(
    df,
    tx_context,
):
    """
    Collapse transcript-specific VEP + GENCODE context
    back to one row per variant.

    Distance filtering itself remains transcript-aware
    through helper pass flags.
    """

    print("\n" + "=" * 60)
    print("COLLAPSE TRANSCRIPT CONTEXT")
    print("=" * 60)

    df = df.copy()

    # General annotation
    df["HGVSc"] = pd.NA
    df["HGVSp"] = pd.NA
    df["Gene_context"] = pd.NA

    # Full transcript-specific distance context
    df["Splice_context"] = pd.NA
    df["TSS_context"] = pd.NA

    # Convenient nearest summaries
    df["Nearest_splice_distance"] = pd.NA
    df["Nearest_splice_transcript"] = pd.NA

    df["Nearest_TSS_distance"] = pd.NA
    df["Nearest_TSS_transcript"] = pd.NA
    df["Nearest_TSS_gene"] = pd.NA

    # Internal filtering helpers
    df["_splice_distance_pass"] = False
    df["_tss_distance_pass"] = False

    grouped = {
        key: group
        for key, group in tx_context.groupby(
            ["Chr_Pos", "Variant_REF_ALT"],
            sort=False,
        )
    }

    for idx, row in df.iterrows():

        key = (
            row["Chr_Pos"],
            row["Variant_REF_ALT"],
        )

        group = grouped.get(key)

        # e.g. mtDNA
        if group is None:
            continue

        hgvsc = []
        hgvsp = []
        gene_context = []
        splice_context = []
        tss_context = []

        # ----------------------------------------------------
        # Pack transcript-specific fields
        # ----------------------------------------------------

        for _, tr in group.iterrows():

            tx = tr["Transcript_ID"]

            if not is_missing(
                tr["HGVSc"]
            ):
                hgvsc.append(
                    f"{tx}={tr['HGVSc']}"
                )

            if not is_missing(
                tr["HGVSp"]
            ):
                hgvsp.append(
                    f"{tx}={tr['HGVSp']}"
                )

            if not is_missing(
                tr["Gene_context"]
            ):
                gene_context.append(
                    f"{tx}="
                    f"{tr['Gene_context']}"
                )

            if pd.notna(
                tr["Splice_distance"]
            ):
                splice_context.append(
                    f"{tx}="
                    f"{int(tr['Splice_distance'])}"
                )

            if pd.notna(
                tr["TSS_distance"]
            ):
                tss_context.append(
                    f"{tx}="
                    f"{int(tr['TSS_distance'])}"
                )

        if hgvsc:
            df.at[idx, "HGVSc"] = (
                " | ".join(
                    dict.fromkeys(hgvsc)
                )
            )

        if hgvsp:
            df.at[idx, "HGVSp"] = (
                " | ".join(
                    dict.fromkeys(hgvsp)
                )
            )

        if gene_context:
            df.at[idx, "Gene_context"] = (
                " | ".join(
                    dict.fromkeys(
                        gene_context
                    )
                )
            )

        if splice_context:
            df.at[idx, "Splice_context"] = (
                " | ".join(
                    splice_context
                )
            )

        if tss_context:
            df.at[idx, "TSS_context"] = (
                " | ".join(
                    tss_context
                )
            )

        # ----------------------------------------------------
        # Nearest splice among SAME Stage 3 transcripts
        # ----------------------------------------------------

        splice = pd.to_numeric(
            group["Splice_distance"],
            errors="coerce",
        )

        if splice.notna().any():

            nearest_idx = (
                splice.abs().idxmin()
            )

            df.at[
                idx,
                "Nearest_splice_distance",
            ] = abs(
                int(
                    splice.loc[
                        nearest_idx
                    ]
                )
            )

            df.at[
                idx,
                "Nearest_splice_transcript",
            ] = group.loc[
                nearest_idx,
                "Transcript_ID",
            ]

        # ----------------------------------------------------
        # Nearest TSS among SAME Stage 3 transcripts
        # ----------------------------------------------------

        tss = pd.to_numeric(
            group["TSS_distance"],
            errors="coerce",
        )

        if tss.notna().any():

            nearest_idx = (
                tss.abs().idxmin()
            )

            df.at[
                idx,
                "Nearest_TSS_distance",
            ] = abs(
                int(
                    tss.loc[
                        nearest_idx
                    ]
                )
            )

            df.at[
                idx,
                "Nearest_TSS_transcript",
            ] = group.loc[
                nearest_idx,
                "Transcript_ID",
            ]

            gene = group.loc[
                nearest_idx,
                "Gene_ID",
            ]

            if not is_missing(gene):
                df.at[
                    idx,
                    "Nearest_TSS_gene",
                ] = gene

        # ----------------------------------------------------
        # Transcript-aware FILTER decisions
        # ----------------------------------------------------

        splice_pass = False

        for _, tr in group.iterrows():

            if (
                "intron_variant"
                not in consequence_terms(
                    tr["Consequence"]
                )
            ):
                continue

            distance = pd.to_numeric(
                tr["Splice_distance"],
                errors="coerce",
            )

            if (
                pd.notna(distance)
                and abs(distance)
                <= SPLICE_DISTANCE_MAX
            ):
                splice_pass = True
                break

        df.at[
            idx,
            "_splice_distance_pass",
        ] = splice_pass

        df.at[
            idx,
            "_tss_distance_pass",
        ] = (
            tss.notna().any()
            and (
                tss.abs()
                <= TSS_DISTANCE_MAX
            ).any()
        )

    print(
        f"Variants annotated: "
        f"{len(df):,}"
    )

    print(
        f"HGVSc: "
        f"{df['HGVSc'].notna().sum():,}"
    )

    print(
        f"HGVSp: "
        f"{df['HGVSp'].notna().sum():,}"
    )

    print(
        f"Gene context: "
        f"{df['Gene_context'].notna().sum():,}"
    )

    print(
        f"Nearest splice: "
        f"{df['Nearest_splice_distance'].notna().sum():,}"
    )

    print(
        f"Nearest TSS: "
        f"{df['Nearest_TSS_distance'].notna().sum():,}"
    )

    print(
        f"Intronic splice pass: "
        f"{df['_splice_distance_pass'].sum():,}"
    )

    print(
        f"TSS <= {TSS_DISTANCE_MAX:,}: "
        f"{df['_tss_distance_pass'].sum():,}"
    )

    return df

def run_distance_filter(df):
    """
    Filter ordinary NEEDS_DISTANCE variants using
    transcript-specific distance decisions already
    calculated in add_distance_annotations().
    """

    splice_pass = (
        df["_intronic"]
        & df["_splice_distance_pass"]
    )

    tss_pass = (
        df["_dels"]
        & df["_tss_distance_pass"]
    )

    regulatory_context = (
        df["Regulatory_context"]
        .fillna("")
        .astype(str)
    )

    nearest_tss = (
        pd.to_numeric(
            df["Nearest_TSS_distance"],
            errors="coerce",
        )
        .abs()
    )

    is_pls = regulatory_context.str.contains(
        r"\bPLS\b",
        case=False,
        regex=True,
    )

    is_pels = regulatory_context.str.contains(
        r"\bpELS\b",
        case=False,
        regex=True,
    )

    pels_relevant = (
        is_pels
        & nearest_tss.le(2000)
    )

    df["Regulatory_relevant"] = (
        is_pls
        | pels_relevant
        | tss_pass
    )

    df["_distance_pass"] = (
        splice_pass
        | tss_pass
    )

    distance_bypass = (
        df["_plp"]
        | df["_case_gene"]
        | df["_second_hit"]
    )

    needs_keep = (
        df["_needs_distance"]
        & (
            distance_bypass
            | df["_distance_pass"]
        )
    )

    keep = (
        df["_direct_keep"]
        | df["_plp"]
        | needs_keep
    )

    print("\n=== Distance filtering ===")

    print(
        f"Input after rarity:                  "
        f"{len(df):,}"
    )

    print(
        f"DIRECT_KEEP (annotation only):       "
        f"{int(df['_direct_keep'].sum()):,}"
    )

    print(
        f"P/LP (annotation only):              "
        f"{int(df['_plp'].sum()):,}"
    )

    print(
        f"Intronic splice <= "
        f"{SPLICE_DISTANCE_MAX} bp:             "
        f"{int(splice_pass.sum()):,}"
    )

    print(
        f"dELS TSS <= "
        f"{TSS_DISTANCE_MAX:,} bp:              "
        f"{int(tss_pass.sum()):,}"
    )

    print(
        f"NEEDS_DISTANCE case/2nd/P-LP bypass: "
        f"{int((df['_needs_distance'] & distance_bypass).sum()):,}"
    )

    print(
        f"NEEDS_DISTANCE passing distance:     "
        f"{int((
            df['_needs_distance']
            & ~distance_bypass
            & df['_distance_pass']
        ).sum()):,}"
    )

    print(
        f"Retained after distance:             "
        f"{int(keep.sum()):,}"
    )

    return df.loc[keep].copy()

# =============================================================================
# Final heterozygous resolution — compound het
# =============================================================================

def resolve_compound_hets(df):
    """
    Mark potential compound-heterozygous candidates after distance filtering.

    Any gene with at least two surviving heterozygous alleles qualifies. The pool
    includes rarity-eligible non-P/LP alleles and heterozygous P/LP alleles, so
    this naturally covers ordinary pairs, P/LP + second hit, and P/LP + P/LP.
    Phase is not proven here; these remain candidate compound-het pairs.
    """
    candidate = df[
        df["_het"]
        & df["_diploid_region"]
        & ~df["_mosaic"]
        & df["VAF"].between(
            GERMLINE_HET_VAF_MIN,
            GERMLINE_HET_VAF_MAX,
            inclusive="both",
        )
        & (df["_het_eligible"] | df["_plp"])
    ]

    # gene -> key -> genomic locus. Different ALT alleles from the same Chr_Pos
    # are not independent hits and therefore cannot form a compound pair.
    by_gene = defaultdict(dict)
    for key, locus, gene_value in zip(
        candidate["_key"], candidate["Chr_Pos"], candidate["Gene_symbol"]
    ):
        for gene in variant_genes(gene_value):
            by_gene[gene][key] = str(locus)

    partners = defaultdict(set)
    genes_by_key = defaultdict(set)
    compound_genes = set()
    for gene, key_to_locus in by_gene.items():
        keys = set(key_to_locus)
        if len(keys) < 2:
            continue
        gene_has_pair = False
        for key in keys:
            valid_partners = {
                partner
                for partner in keys - {key}
                if key_to_locus[partner] != key_to_locus[key]
            }
            if not valid_partners:
                continue
            gene_has_pair = True
            genes_by_key[key].add(gene)
            partners[key].update(valid_partners)
        if gene_has_pair:
            compound_genes.add(gene)

    df["_compound_het"] = df["_key"].isin(partners)
    df["Compound_het_partner"] = df["_key"].map(
        lambda key: TX_SEP.join(sorted(partners.get(key, ())))
    )
    df["Compound_het_gene"] = df["_key"].map(
        lambda key: TX_SEP.join(sorted(genes_by_key.get(key, ())))
    )

    print("\n=== Compound heterozygous ===")
    print(f"Compound candidate variants: {int(df['_compound_het'].sum()):,}")
    print(f"Genes with compound pairs:   {len(compound_genes):,}")
    return df

# =============================================================================
# Final heterozygous resolution — singleton AD
# =============================================================================

def run_final_heterozygous_filter(df):
    """
    Heterozygous candidates without a compound partner must pass
    the stricter AD singleton thresholds.

    Ordinary gene:
        AF <= 1e-4
        gnomAD het < 30
        gnomAD hom = 0

    Case gene:
        AF <= 1e-3
        gnomAD het < 50
        gnomAD hom = 0

    P/LP remains retained independently.
    """

    candidate = (
        df["_het_eligible"]
        & df["_het"]
        & ~df["_compound_het"]
    )

    ordinary_rule = RARITY_THRESHOLDS["AD_singleton"]
    case_rule = RARITY_THRESHOLDS["AD_case_singleton"]

    # ---------------------------------------------------------------------
    # Ordinary singleton
    # ---------------------------------------------------------------------

    ordinary = (
        candidate
        & ~df["_case_gene"]
        & (
            df["_filter_af"].le(ordinary_rule["af_max"])
            | df["_filter_af"].isna()
        )
        & (
            df["gnomAD_het"].lt(ordinary_rule["het_lt"])
            | df["gnomAD_het"].isna()
        )
        & (
            df["gnomAD_hom"].le(ordinary_rule["hom_max"])
            | df["gnomAD_hom"].isna()
        )
    )

    # ---------------------------------------------------------------------
    # Case-gene singleton
    # ---------------------------------------------------------------------

    case = (
        candidate
        & df["_case_gene"]
        & (
            df["_filter_af"].le(case_rule["af_max"])
            | df["_filter_af"].isna()
        )
        & (
            df["gnomAD_het"].lt(case_rule["het_lt"])
            | df["gnomAD_het"].isna()
        )
        & (
            df["gnomAD_hom"].le(case_rule["hom_max"])
            | df["gnomAD_hom"].isna()
        )
    )

    df["_ad_singleton"] = ordinary | case

    print("\n=== Final heterozygous filtering ===")
    print(f"Compound variants:       {int(df['_compound_het'].sum()):,}")
    print(f"Ordinary AD singletons:  {int(ordinary.sum()):,}")
    print(f"Case-gene AD singletons: {int(case.sum()):,}")
    print(f"P/LP retained:           {int(df['_plp'].sum()):,}")
    print(f"Secondary retained:      {int(df['_secondary_reportable'].sum()):,}")

    return df


def inheritance_categories(row, sex):
    """Return all supported inheritance routes for one retained variant."""
    categories = []

    if row.get("_compound_het", False):
        categories.append("AR_compound_heterozygous")
    if row.get("_ad_singleton", False):
        categories.append(
            "X_linked_heterozygous"
            if row["_chrom"] == "X" and sex == "female" and not row["_is_par"]
            else "AD"
        )
    if row.get("_hom", False):
        categories.append(
            "X_linked_homozygous"
            if row["_chrom"] == "X" and sex == "female" and not row["_is_par"]
            else "AR_homozygous"
        )
    if row.get("_x_hemi", False):
        categories.append("X_linked_hemizygous")
    if row.get("_mosaic", False):
        categories.append("Mosaic")
    if row.get("_mito", False):
        categories.append("Mitochondrial")

    # A P/LP variant bypasses rarity, so if no route above resolved its
    # inheritance, infer only from an observed ALT genotype; never label an
    # unknown/partial/no-ALT call as AD by default.
    if row.get("_plp", False) and not categories:
        if row["_chrom"] == "MT" and row["_has_alt"]:
            categories.append("Mitochondrial")
        elif row["_chrom"] == "X" and sex == "male" and not row["_is_par"] and (
            row["_hom_alt"] or row["_haploid_alt"]
        ):
            categories.append("X_linked_hemizygous")
        elif row["_chrom"] == "X" and sex == "female" and not row["_is_par"]:
            if row["_hom_alt"]:
                categories.append("X_linked_homozygous")
            elif row["_het"]:
                categories.append("X_linked_heterozygous")
        elif row["_hom_alt"]:
            categories.append("AR_homozygous")
        elif row["_het"]:
            categories.append("AD")

    if row.get("_secondary_reportable", False):
        categories.append("Secondary")

    return TX_SEP.join(dict.fromkeys(categories))


def build_final_output(df, sex):
    """
    Keep every variant surviving at least one final Stage-3 route
    and return the clean user-facing table.
    """

    final_keep = (
        df["_hom"]
        | df["_x_hemi"]
        | df["_mosaic"]
        | df["_mito"]
        | df["_compound_het"]
        | df["_ad_singleton"]
        | df["_plp"]
    )

    out = df.loc[final_keep].copy()

    out["Case_gene_match"] = out["_case_gene"]
    out["Secondary_gene_match"] = out["_secondary_gene"]
    out["Secondary_MOI"] = out["_secondary_moi"]
    out["Secondary_reportable"] = out["_secondary_reportable"]

    # ------------------------------------------------------------------
    # Inheritance category
    # ------------------------------------------------------------------

    out["inheritance_category"] = out.apply(
        inheritance_categories,
        axis=1,
        sex=sex,
    )

    print("\n=== Final Stage 3 ===")
    print(f"Input after distance: {len(df):,}")
    print(f"Final variants:       {len(out):,}")

    return out[FINAL_COLUMNS]

# =============================================================================
# Validation and output
# =============================================================================

def validate_final(final):
    """Fail fast on structural, inheritance, Secondary, and compound-pair errors."""
    validate_table(final, "Stage 3 final", FINAL_COLUMNS, aligned=True)
    if final.empty:
        raise ValueError("Stage 3 produced zero final variants.")

    missing_cols = [c for c in FINAL_COLUMNS if c not in final.columns]
    if missing_cols:
        raise ValueError(f"Final table is missing columns: {missing_cols}")

    if final[KEY_COLS].duplicated().any():
        dup = final.loc[final[KEY_COLS].duplicated(keep=False), KEY_COLS].head()
        raise ValueError(f"Duplicate variants in final Stage-3 table:\n{dup}")

    categories = final["inheritance_category"].fillna("").astype(str)
    if categories.str.strip().eq("").any():
        n = int(categories.str.strip().eq("").sum())
        raise ValueError(f"{n} final variants have no inheritance_category.")

    category_sets = categories.map(
        lambda value: {x.strip() for x in value.split(TX_SEP) if x.strip()}
    )
    secondary_label = category_sets.map(lambda x: "Secondary" in x)
    secondary_flag = final["Secondary_reportable"].fillna(False).astype(bool)
    if not secondary_label.equals(secondary_flag):
        raise ValueError("Secondary_reportable and inheritance_category disagree.")

    keys = final["Chr_Pos"].astype(str) + "|" + final["Variant_REF_ALT"].astype(str)
    key_set = set(keys)
    locus_by_key = dict(zip(keys, final["Chr_Pos"].astype(str)))
    partner_map = {}
    compound_rows = category_sets.map(lambda x: "AR_compound_heterozygous" in x)

    for key, partner_value, gene_value, is_compound in zip(
        keys, final["Compound_het_partner"], final["Compound_het_gene"], compound_rows
    ):
        partners = {
            x for x in split_tx_field(partner_value)
            if not is_missing(x)
        }
        partner_map[key] = partners

        if is_compound and (not partners or is_missing(gene_value)):
            raise ValueError(f"Compound-het variant lacks partner/gene: {key}")
        if key in partners:
            raise ValueError(f"Compound-het variant lists itself as partner: {key}")
        missing_partners = partners - key_set
        if missing_partners:
            raise ValueError(f"Compound partner(s) absent from final table for {key}: {missing_partners}")
        same_locus_partners = {
            partner for partner in partners
            if locus_by_key.get(partner) == locus_by_key.get(key)
        }
        if same_locus_partners:
            raise ValueError(
                f"Compound partner(s) share the same genomic locus for {key}: "
                f"{same_locus_partners}"
            )

    for key, partners in partner_map.items():
        for partner in partners:
            if key not in partner_map.get(partner, set()):
                raise ValueError(f"Non-reciprocal compound pair: {key} <-> {partner}")

    print("\n=== Final validation ===")
    print(f"Unique final variants: {len(final):,}")
    print("Validation: PASS")


def write_filtered_keys(final, path):
    """
    Write only the variant keys retained by Stage 3.
    """

    keys = (
        final[KEY_COLS]
        .drop_duplicates()
        .copy()
    )

    keys.to_csv(
        path,
        sep="\t",
        index=False,
        compression="gzip",
    )

    print(f"Filtered keys written: {len(keys):,}")


def write_final_csv(final, path):
    final.to_csv(
        path,
        index=False,
    )

    print(f"Final CSV written:     {len(final):,}")


def print_final_summary(stage_counts, final):
    """
    Print:
      - how many variants survived each stage
      - how many were removed at each stage
      - counts by inheritance category
    """

    print("\n" + "=" * 60)
    print("STAGE 3 FINAL SUMMARY")
    print("=" * 60)

    print("\nFiltering flow")
    print("-" * 60)

    previous = None

    for stage, count in stage_counts.items():

        if previous is None:
            removed = 0
        else:
            removed = previous - count

        print(
            f"{stage:<34}"
            f"{count:>12,}"
            f"   removed: {removed:>10,}"
        )

        previous = count

    print("\nInheritance categories")
    print("-" * 60)

    category_counts = (
        final["inheritance_category"]
        .fillna("")
        .astype(str)
        .str.split(TX_SEP, regex=False)
        .explode()
        .str.strip()
    )

    category_counts = category_counts[
        category_counts.ne("")
    ].value_counts()

    for category, count in category_counts.items():
        print(
            f"{category:<34}"
            f"{count:>12,}"
        )

    print("\nFinal unique variants")
    print("-" * 60)
    print(f"{len(final):,}")

def main():
    started = time.perf_counter()
    stage_counts = {}

    # ---------------------------------------------------------------------
    # 1. Load + QC
    # ---------------------------------------------------------------------
    df, n_collapsed = load_collapsed(
        str(snakemake.input.collapsed)
    )

    stage_counts["Collapsed input"] = n_collapsed
    stage_counts["After QC"] = len(df)
    input_keys = set(variant_keys(df, "Stage 3 input"))

    # 2. Functional classification
    # GT flags are needed by later inheritance and Secondary logic.
    sex, case_genes, hom_max_override = load_case_context(str(snakemake.input.case_context))
    df = add_genotype_flags(df)
    df = run_functional_classification(df)
    stage_counts["After functional filtering"] = len(df)

    # 3. Case / P-LP / ACMG Secondary context
    df = add_case_and_genotype_context(df, sex, case_genes)
    acmg_moi = load_acmg_sf(str(snakemake.input.acmg_sf))
    df = add_secondary_context(df, sex, acmg_moi)
    # 4. Initial genotype + rarity
    df = run_initial_rarity_filters(df, sex, hom_max_override)
    stage_counts["After genotype + rarity"] = len(df)

    # 5. Distance annotation + filtering
    vep_results, id_to_variant = (
        load_or_run_stage3_vep(
            df,
            STAGE3_VEP_CACHE,
        )
    )
    tx_context = build_transcript_context(
        df,
        vep_results,
        id_to_variant,
        snakemake.input.gencode_gtf,
    )

    df = add_distance_annotations(
        df,
        tx_context,
    )

    df = run_distance_filter(df)
    stage_counts["After distance filtering"] = len(df)

    print("\n=== After distance filtering ===")
    print(f"Active variants: {len(df):,}")

    # 6. Final heterozygous resolution
    df = resolve_compound_hets(df)
    df = run_final_heterozygous_filter(df)

    # 7. Final routes + clean output table
    final = build_final_output(df, sex)
    stage_counts["Final retained"] = len(final)
    # 8. Validation and outputs
    validate_final(final)
    if not set(variant_keys(final, "Stage 3 output")).issubset(input_keys):
        raise ValueError("Stage 3 added unexpected variants")
    record_qc("stage3_filter", n_collapsed, len(final), n_collapsed, len(final), "filter", counts=stage_counts)
    write_final_csv(
        final,
        str(snakemake.output.candidates),
    )
    write_filtered_keys(
        final,
        str(snakemake.output.filtered_keys),
    )

    print_final_summary(stage_counts, final)
    elapsed = time.perf_counter() - started
    print(f"\nStage 3 runtime: {elapsed / 60:.2f} min ({elapsed:.1f} s)")

if __name__ == "__main__":
    main()
