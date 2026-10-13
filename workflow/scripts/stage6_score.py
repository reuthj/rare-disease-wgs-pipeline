#!/usr/bin/env python3
"""Stage 6: weighted candidate scoring and prioritization."""

import argparse
import math
import os
import re
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, NamedTuple

from main_classification import is_main_plp, classification_tokens
import pandas as pd
from workflow_checks import validate_table, record_qc, is_mitochondrial_variant
import yaml


STAGE6_OUTPUT_COLUMNS = [
    "Variant_pathogenicity_score",
    "Pathogenicity_evidence_coverage",
    "Inheritance_score",
    "Gene_disease_evidence_score",
    "Case_relevance_score",
    "Rarity_score",
    "Scored_gene",
    "Scored_transcript",
    "Scored_route",
    "Scored_inheritance",
    "AD_fallback",
    "Scored_compound_partner",
    "Compound_pair_score",
    "Final_score",
]


REQUIRED_STAGE5_COLUMNS = [
    "Chr_Pos",
    "Variant_REF_ALT",
    "Gene_symbol",
    "Gene_symbol_transcript_aligned",
    "Transcript_ID",
    "inheritance_category",
    "Case_gene_match",
    "ClinVar_significance",
    "CADD",
    "phyloP100way",
    "phastCons100way",
    "Disease_associations",
    "Missense_constraint",
    "LOEUF",
    "ClinVar_pathogenic_fraction",
    "ClinVar_pathogenic_count",
    "Known_variant_landscape",
    "Gene_HPO_terms",
    "Matched_HPO_terms",
    "Inheritance_and_mechanism",
    "Literature_support",
    "Top_publications",
    "MANE_Transcript_ID",
]


REQUIRED_SCORING_COLUMNS = [
    "Transcript_biotype",
    "Consequences",
    "HGVSc",
    "HGVSp",
    "Regulatory_context",
    "Splice_context",
    "Nearest_splice_distance",
    "Nearest_splice_transcript",
    "TSS_context",
    "Nearest_TSS_distance",
    "Nearest_TSS_transcript",
    "Nearest_TSS_gene",
    "Regulatory_relevant",
    "Compound_het_gene",
    "Compound_het_partner",
    "gnomAD_AF",
    "gnomAD_grpmax_AF",
    "gnomAD_het",
    "gnomAD_hom",
    "gnomAD_hemi",
    "gnomAD_mt_AF_hom",
    "gnomAD_mt_AF_het",
    "is_missense",
    "is_lof",
    "is_inframe_indel",
    "is_splice_relevant",
    "is_regulatory",
    "is_mitochondrial",
    "Evidence_routes",
    "Transcript_routes",
    "Protein_context",
    "NMD_context",
    "Missense_evidence",
    "Splice_prediction",
    "Regulatory_prediction",
    "LoF_evidence",
    "Inframe_indel_evidence",
    "Mitochondrial_evidence",
]


GENE_EVIDENCE_FIELDS = [
    "Disease_associations",
    "Missense_constraint",
    "LOEUF",
    "ClinVar_pathogenic_fraction",
    "ClinVar_pathogenic_count",
    "Known_variant_landscape",
    "Gene_HPO_terms",
    "Matched_HPO_terms",
    "Inheritance_and_mechanism",
    "Literature_support",
    "Top_publications",
]


# Only evidence that is expected to be reusable for the same gene is cached.
# Literature fields stay row-specific because they may include variant-specific
# publications even when the text also contains a gene name.
GENE_CACHE_TEXT_FIELDS = [
    "Disease_associations",
    "Missense_constraint",
    "Known_variant_landscape",
    "Gene_HPO_terms",
    "Matched_HPO_terms",
    "Inheritance_and_mechanism",
]


GENE_CACHE_NUMERIC_FIELDS = [
    "LOEUF",
    "ClinVar_pathogenic_fraction",
    "ClinVar_pathogenic_count",
]


GENE_CACHE_FIELDS = GENE_CACHE_TEXT_FIELDS + GENE_CACHE_NUMERIC_FIELDS


TRANSCRIPT_EVIDENCE_FIELDS = [
    "HGVSc",
    "HGVSp",
    "Splice_context",
    "TSS_context",
    "Protein_context",
    "NMD_context",
    "Missense_evidence",
    "Splice_prediction",
    "LoF_evidence",
    "Inframe_indel_evidence",
]


REQUIRED_CASE_FIELDS = [
    "sex",
    "expected_inheritance",
    "disease_names",
    "hpo_terms",
    "broad_disease_domains",
    "manual_genes",
    "case_genes",
]


LIST_CASE_FIELDS = [
    "expected_inheritance",
    "disease_names",
    "hpo_terms",
    "broad_disease_domains",
    "manual_genes",
    "case_genes",
]


EXPECTED_INHERITANCE_MODELS = {
    "AD",
    "AD_de_novo",
    "AR_homozygous",
    "AR_compound_heterozygous",
    "X_linked",
    "Mosaic",
    "Mitochondrial",
}


INHERITANCE_ALIASES = {
    "AD_singleton": "AD",
    "X_linked_hemizygous": "X_linked",
    "X_linked_hemi": "X_linked",
    "mtDNA": "Mitochondrial",
}


NO_GENE_SCORED_GENE = ""
NO_GENE_GENE_DISEASE_SCORE = 0.0
NO_GENE_CASE_RELEVANCE_SCORE = 0.0


class Stage6InputError(ValueError):
    """Raised when Stage 6 inputs are structurally unsafe to score."""


MISSING_TOKENS = {"", "na", "n/a", "nan", "none", "null", ".", "-"}
ALIGNED_SEPARATOR = re.compile(r"\s*\|\s*")
NUMBER_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"


CADD_ANCHORS = [(5.0, 0.0), (10.0, 2.0), (15.0, 4.0), (20.0, 6.0), (22.5, 8.0), (25.0, 10.0)]
PHYLOP_ANCHORS = [(0.0, 0.0), (0.5, 2.0), (1.0, 4.0), (1.5, 6.0), (2.0, 8.0), (2.5, 10.0)]
PHASTCONS_ANCHORS = [(0.0, 0.0), (0.15, 2.0), (0.30, 4.0), (0.45, 6.0), (0.60, 8.0), (0.75, 10.0)]
REVEL_ANCHORS = [(0.18, 0.0), (0.30, 2.0), (0.50, 4.0), (0.64, 6.0), (0.77, 8.0), (0.93, 10.0)]
ALPHAMISSENSE_ANCHORS = [(0.15, 0.0), (0.30, 4.0), (0.45, 6.0), (0.60, 8.0), (0.75, 10.0)]
SPLICEAI_ANCHORS = [(0.0, 0.0), (0.10, 2.0), (0.20, 4.0), (0.50, 6.0), (0.65, 8.0), (0.80, 10.0)]


RARITY_AF_THRESHOLDS = {
    "AD": 1e-4,
    "AD_case_gene": 1e-3,
    "AR_compound_heterozygous": 1e-3,
    "AR_compound_heterozygous_case_gene": 1e-2,
    "AR_homozygous": 1e-3,
    "X_linked": 5e-3,
    "Mosaic": 1e-3,
    "Mitochondrial": 5e-3,
}


RARITY_LOG_RATIO_ANCHORS = [
    (-1.0, 0.0),
    (0.0, 1.0),
    (math.log10(5.0), 5.0),
    (math.log10(50.0), 8.0),
    (math.log10(500.0), 10.0),
]


class TranscriptPathogenicity(NamedTuple):
    transcript: str
    gene: str | None
    route: str
    score: float
    coverage: float


class Candidate(NamedTuple):
    """One complete scoring hypothesis for one variant and one gene."""

    gene: str | None
    gene_order: int
    transcript: str
    pathogenicity_route: str
    variant_pathogenicity_score: float
    pathogenicity_evidence_coverage: float
    inheritance_score: float
    gene_disease_evidence_score: float
    case_relevance_score: float
    rarity_score: float
    final_score: float
    is_case_gene: bool
    is_pathogenic_clinvar: bool
    scored_compound_partner: str = ""
    compound_pair_score: float | None = None
    inheritance_model: str = ""
    ad_fallback: bool = False

    @property
    def has_priority_exception(self) -> bool:
        return self.is_case_gene or self.is_pathogenic_clinvar


class CompoundPair(NamedTuple):
    """A Stage 3-defined compound pair scored as one ranking unit."""

    gene: str
    left_index: int
    right_index: int
    left_candidate: Candidate
    right_candidate: Candidate
    pair_score: float


def is_missing(value: Any) -> bool:
    if value is None:
        return True
    try:
        if pd.isna(value):
            return True
    except (TypeError, ValueError):
        pass
    return str(value).strip().lower() in MISSING_TOKENS


def to_float(value: Any) -> float | None:
    if is_missing(value):
        return None
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return None


def linear_score(
    value: float | None,
    anchors: Sequence[tuple[float, float]],
    missing_score: float = 5.0,
) -> float:
    """Linearly interpolate a 0-10 evidence score between locked anchors."""
    if value is None:
        return missing_score
    if not anchors:
        raise ValueError("linear_score requires at least one anchor")

    ordered = sorted(anchors)
    if value <= ordered[0][0]:
        return ordered[0][1]
    if value >= ordered[-1][0]:
        return ordered[-1][1]

    for (left_x, left_y), (right_x, right_y) in zip(
        ordered,
        ordered[1:],
        strict=False,
    ):
        if left_x <= value <= right_x:
            fraction = (value - left_x) / (right_x - left_x)
            return left_y + fraction * (right_y - left_y)

    raise RuntimeError("Could not interpolate value between anchors")


def extract_labeled_float(value: Any, label: str) -> float | None:
    """Extract `label=number` from one Stage 5 evidence string."""
    if is_missing(value):
        return None
    match = re.search(
        rf"(?:^|[;|])\s*{re.escape(label)}\s*=\s*({NUMBER_PATTERN})",
        str(value),
        flags=re.IGNORECASE,
    )
    return to_float(match.group(1)) if match else None


def parse_cadd_phred(value: Any) -> float | None:
    """Read the PHRED value from Stage 5 CADD text, or accept a bare number."""
    labeled = extract_labeled_float(value, "PHRED")
    if labeled is not None:
        return labeled
    return to_float(value)


def score_cadd(value: Any) -> tuple[float, bool]:
    """Return CADD points out of 20 and whether CADD was observed."""
    parsed = parse_cadd_phred(value)
    return linear_score(parsed, CADD_ANCHORS) * 2.0, parsed is not None


def score_conservation(
    phylop_value: Any,
    phastcons_value: Any,
) -> tuple[float, int]:
    """Return conservation points out of 10 and number of observed sources."""
    phylop = to_float(phylop_value)
    phastcons = to_float(phastcons_value)
    observed_scores: list[float] = []

    if phylop is not None:
        observed_scores.append(linear_score(phylop, PHYLOP_ANCHORS))
    if phastcons is not None:
        observed_scores.append(linear_score(phastcons, PHASTCONS_ANCHORS))

    if not observed_scores:
        return 5.0, 0
    return sum(observed_scores) / len(observed_scores), len(observed_scores)


def to_bool(value: Any) -> bool | None:
    if is_missing(value):
        return None
    normalized = str(value).strip().casefold()
    if normalized in {"true", "yes", "1"}:
        return True
    if normalized in {"false", "no", "0"}:
        return False
    return None


def score_splice_distance(value: Any) -> tuple[float, bool]:
    distance = to_float(value)
    if distance is None:
        return 5.0, False
    distance = abs(distance)
    if distance <= 5:
        return 10.0, True
    if distance <= 20:
        return 5.0, True
    return 0.0, True


def score_loeuf(value: float | None) -> tuple[float, bool]:
    if value is None:
        return 5.0, False
    if value < 0.35:
        return 10.0, True
    if value <= 0.90:
        return 5.0, True
    return 0.0, True


def score_tss_distance(value: Any) -> tuple[float, bool]:
    distance = to_float(value)
    if distance is None:
        return 5.0, False
    distance = abs(distance)
    if distance <= 1_000:
        return 10.0, True
    if distance <= 2_000:
        return 8.0, True
    if distance <= 5_000:
        return 6.0, True
    if distance <= 10_000:
        return 4.0, True
    return 2.0, True


def score_remm(value: Any) -> tuple[float, bool]:
    remm = extract_labeled_float(value, "ReMM")
    if remm is None:
        remm = to_float(value)
    if remm is None:
        return 5.0, False
    if remm < 0.50:
        return 0.0, True
    if remm < 0.75:
        return 2.0, True
    if remm < 0.90:
        return 4.0, True
    if remm < 0.95:
        return 6.0, True
    if remm < 0.98:
        return 8.0, True
    return 10.0, True


def score_protein_context(value: Any) -> tuple[float, bool]:
    if is_missing(value):
        return 5.0, False
    text = str(value).strip()
    domain_match = re.search(r"(?:^|;)\s*Domains\s*=\s*([^;]+)", text, re.I)
    if domain_match and not is_missing(domain_match.group(1)):
        return 10.0, True
    if re.search(r"(?:^|;)\s*(?:AA|Protein_pos)\s*=", text, re.I):
        return 5.0, True
    return 0.0, True


def score_lof_context(lof_value: Any, nmd_value: Any) -> tuple[float, bool]:
    if is_missing(lof_value) and is_missing(nmd_value):
        return 5.0, False
    text = f"{lof_value};{nmd_value}".casefold()
    if "loftee=hc" in text or re.search(r"nmd[_ -]?(?:predicted|trigger)", text):
        return 10.0, True
    if "loftee=lc" in text or "nmd_escaping" in text or "nmd escaping" in text:
        return 5.0, True
    if any(term in text for term in ("not_lof", "not lof", "no_lof", "no lof")):
        return 0.0, True
    return 5.0, True


def score_inframe_evidence(value: Any) -> tuple[float, bool]:
    if is_missing(value):
        return 5.0, False
    text = str(value).casefold()
    clear_terms = ("active_site", "binding_site", "functional_site", "deleterious")
    absent_terms = ("no_support", "not_supporting", "benign")
    if any(term in text for term in clear_terms):
        return 10.0, True
    if any(term in text for term in absent_terms):
        return 0.0, True
    return 5.0, True


def score_regulatory_annotation(
    regulatory_context: Any,
    regulatory_relevant: Any,
) -> tuple[float, bool]:
    relevant = to_bool(regulatory_relevant)
    if relevant is False:
        return 0.0, True
    if is_missing(regulatory_context):
        return 5.0, False
    text = str(regulatory_context).casefold()
    if "promoter" in text or re.search(r"(?:^|:)pls(?:\b|$)", text):
        return 10.0, True
    if any(term in text for term in ("pels", "dels", "enhancer")):
        return 5.0, True
    return 5.0 if relevant is True else 0.0, True


def clinvar_pathogenicity_adjustment(value: Any) -> float:
    if is_missing(value):
        return 0.0
    tokens_main = classification_tokens(value)
    if is_main_plp(value):
        return 20.0 if all(t in {"p", "pathogenic"} for t in tokens_main) else 15.0
    if any(t in {"p", "lp", "pathogenic", "likely pathogenic"} for t in tokens_main):
        return 0.0
    normalized = str(value).casefold().replace("_", " ")
    tokens = [
        token.strip()
        for token in re.split(r"[|/;,]", normalized)
        if token.strip()
    ]
    has_conflicting = any("conflict" in token for token in tokens)
    has_other = any(token == "other" for token in tokens)
    has_pathogenic = any(token == "pathogenic" for token in tokens)
    has_likely_pathogenic = any(token == "likely pathogenic" for token in tokens)
    has_benign = any(token == "benign" for token in tokens)
    has_likely_benign = any(token == "likely benign" for token in tokens)
    has_vus = any(token == "uncertain significance" for token in tokens)
    if has_pathogenic:
        if has_conflicting or has_other or has_vus or has_benign or has_likely_benign:
            return 10.0
        if has_likely_pathogenic:
            return 15.0
        return 20.0
    if has_likely_pathogenic:
        return 10.0 if (
            has_conflicting or has_other or has_vus or has_benign or has_likely_benign
        ) else 15.0
    if has_benign:
        if has_conflicting or has_other or has_vus:
            return -10.0
        if has_likely_benign:
            return -15.0
        return -20.0
    if has_likely_benign:
        return -10.0 if (has_conflicting or has_other or has_vus) else -15.0
    if has_vus:
        return 0.0 if (has_conflicting or has_other) else 5.0
    return 0.0

def clinvar_has_pathogenic_evidence(value: Any) -> bool:
    """Use the authoritative main classification only."""
    return is_main_plp(value)

def clamp_score(value: float) -> float:
    return min(100.0, max(0.0, value))


def split_aligned(value: Any) -> list[str]:
    """Split a positional ` | ` field without removing empty positions."""
    if value is None:
        return []
    text = str(value).strip().replace(r"\|", "|")
    if not text:
        return []
    return [part.strip() for part in ALIGNED_SEPARATOR.split(text)]


def _stable_unique(values: Sequence[str]) -> list[str]:
    return list(dict.fromkeys(values))


def parse_gene_symbols(value: Any) -> list[str]:
    genes = [
        gene
        for gene in split_aligned(value)
        if not is_missing(gene)
    ]
    return _stable_unique(genes)


def parse_row_genes(row: Mapping[str, Any]) -> list[str]:
    return parse_gene_symbols(row.get("Gene_symbol"))


def parse_transcript_gene_map(
    row: Mapping[str, Any],
) -> dict[str, str | None]:
    transcripts = split_aligned(row.get("Transcript_ID"))
    aligned_genes = split_aligned(row.get("Gene_symbol_transcript_aligned"))

    if not transcripts:
        if is_mitochondrial_variant(row):
            genes = parse_row_genes(row)
            if len(aligned_genes) > 1:
                raise Stage6InputError("MT variant without transcripts has contradictory aligned genes.")
            aligned_gene = aligned_genes[0] if aligned_genes else None
            if is_missing(aligned_gene):
                aligned_gene = genes[0] if len(genes) == 1 else None
            elif aligned_gene not in genes:
                raise Stage6InputError("MT aligned gene is absent from Gene_symbol.")
            # An empty transcript remains empty in the output. COMMON evidence
            # supports a variant-level hypothesis without inventing a transcript.
            return {"": aligned_gene}
        raise Stage6InputError("A row has no Transcript_ID values.")
    if len(transcripts) != len(aligned_genes):
        key = f"{row.get('Chr_Pos', '?')}|{row.get('Variant_REF_ALT', '?')}"
        raise Stage6InputError(
            "Transcript/gene alignment length mismatch for "
            f"{key}: {len(transcripts)} transcripts and "
            f"{len(aligned_genes)} aligned genes."
        )
    if any(is_missing(transcript) for transcript in transcripts):
        raise Stage6InputError("Transcript_ID contains an empty aligned position.")
    if len(set(transcripts)) != len(transcripts):
        raise Stage6InputError("Transcript_ID contains duplicates within one row.")

    row_genes = parse_row_genes(row)
    canonical_by_casefold: dict[str, str] = {}
    for row_gene in row_genes:
        key = row_gene.casefold()
        if key in canonical_by_casefold and canonical_by_casefold[key] != row_gene:
            raise Stage6InputError(
                "Gene_symbol contains symbols that differ only by letter case: "
                f"{canonical_by_casefold[key]!r}, {row_gene!r}."
            )
        canonical_by_casefold[key] = row_gene

    mapping: dict[str, str | None] = {}
    for transcript, gene in zip(transcripts, aligned_genes, strict=True):
        clean_gene = None
        if not is_missing(gene):
            clean_gene = canonical_by_casefold.get(gene.casefold())
            if clean_gene is None:
                raise Stage6InputError(
                    f"Aligned gene {gene!r} is absent from Gene_symbol."
                )
        mapping[transcript] = clean_gene
    return mapping


def parse_gene_blocks(
    value: Any,
    genes: Sequence[str],
    warnings: list[str] | None = None,
    field_name: str = "gene evidence",
) -> dict[str, str | None]:
    """Parse gene-prefixed blocks without splitting internal `|` or `||`."""
    unique_genes = _stable_unique([gene for gene in genes if not is_missing(gene)])
    result = {gene: None for gene in unique_genes}
    if not unique_genes or is_missing(value):
        return result

    text = str(value).strip().replace(r"\|", "|").replace("\u00a0", " ")
    gene_pattern = "|".join(
        re.escape(gene) for gene in sorted(unique_genes, key=len, reverse=True)
    )
    marker = re.compile(
        rf"(?P<prefix>^|\s*\|\|\s*)(?P<gene>{gene_pattern})\s*:\s*"
    )
    matches = list(marker.finditer(text))

    if not matches:
        if len(unique_genes) == 1:
            result[unique_genes[0]] = text
        elif warnings is not None:
            warnings.append(
                f"{field_name}: unprefixed evidence for a multi-gene row"
            )
        return result

    leading_text = text[: matches[0].start()].strip()
    if leading_text and warnings is not None:
        warnings.append(f"{field_name}: ignored text before first gene block")

    for index, match in enumerate(matches):
        gene = match.group("gene")
        end = matches[index + 1].start() if index + 1 < len(matches) else len(text)
        evidence = text[match.end() : end].strip()
        evidence = re.sub(r"\s*\|\|\s*$", "", evidence).strip()
        parsed = None if is_missing(evidence) else evidence

        if result[gene] is not None and parsed is not None:
            if warnings is not None:
                warnings.append(f"{field_name}: repeated block for {gene}")
            result[gene] = f"{result[gene]} || {parsed}"
        elif parsed is not None:
            result[gene] = parsed

    return result


def parse_tx_blocks(
    value: Any,
    transcripts: Sequence[str] | None = None,
    warnings: list[str] | None = None,
    field_name: str = "transcript evidence",
) -> dict[str, str | None]:
    """Parse `ENST...=value | ENST...=value` evidence blocks."""
    result: dict[str, str | None] = {}
    if is_missing(value):
        return result

    for block in split_aligned(value):
        if is_missing(block):
            continue
        if "=" not in block:
            if warnings is not None:
                warnings.append(f"{field_name}: block without '=': {block}")
            continue
        raw_transcript, evidence = block.split("=", 1)
        raw_transcript = raw_transcript.strip()
        if not raw_transcript:
            if warnings is not None:
                warnings.append(f"{field_name}: empty transcript identifier")
            continue
        transcript = raw_transcript
        if transcripts is not None:
            transcript = resolve_transcript_id(raw_transcript, transcripts)
            if transcript is None:
                raise Stage6InputError(
                    f"{field_name}: evidence references unknown transcript "
                    f"{raw_transcript!r}."
                )
        parsed = None if is_missing(evidence) else evidence.strip()
        if transcript in result and warnings is not None:
            warnings.append(f"{field_name}: repeated block for {transcript}")
        result[transcript] = parsed
    return result


def transcript_base_id(transcript: str) -> str:
    return transcript.strip().split(".", 1)[0]


def resolve_transcript_id(
    identifier: str,
    transcripts: Sequence[str],
) -> str | None:
    clean_identifier = identifier.strip()
    if clean_identifier in transcripts:
        return clean_identifier

    base_identifier = transcript_base_id(clean_identifier)
    matches = [
        transcript
        for transcript in transcripts
        if transcript_base_id(transcript) == base_identifier
    ]
    if len(matches) == 1:
        return matches[0]
    return None


def parse_routes(value: Any) -> list[str]:
    if is_missing(value):
        return []
    routes = re.split(r"\s*(?:\||,)\s*", str(value).strip())
    return _stable_unique([route.upper() for route in routes if not is_missing(route)])


def parse_transcript_routes(
    value: Any,
    transcripts: Sequence[str],
) -> dict[str, list[str]]:
    result = {transcript: [] for transcript in transcripts}
    if is_missing(value):
        return result

    for block in split_aligned(value):
        if is_missing(block):
            continue
        if ":" not in block:
            raise Stage6InputError(
                f"Transcript_routes block has no ':' separator: {block!r}."
            )
        raw_transcript, raw_routes = block.split(":", 1)
        transcript = resolve_transcript_id(raw_transcript, transcripts)
        if transcript is None:
            raise Stage6InputError(
                "Transcript_routes references unknown transcript "
                f"{raw_transcript.strip()!r}."
            )
        if result[transcript]:
            raise Stage6InputError(
                f"Transcript_routes repeats transcript {transcript!r}."
            )
        result[transcript] = parse_routes(raw_routes)
    return result


def parse_mane_transcripts(
    value: Any,
    transcripts: Sequence[str],
) -> list[str]:
    if is_missing(value):
        return []

    mane_transcripts: list[str] = []
    for raw_transcript in split_aligned(value):
        if is_missing(raw_transcript):
            continue
        transcript = resolve_transcript_id(raw_transcript, transcripts)
        if transcript is None:
            raise Stage6InputError(
                "MANE_Transcript_ID references unknown transcript "
                f"{raw_transcript!r}."
            )
        mane_transcripts.append(transcript)
    return _stable_unique(mane_transcripts)


def parse_inheritance_categories(value: Any) -> list[str]:
    if is_missing(value):
        return []
    parts = re.split(r"\s*(?:\|\||\||,|;)\s*", str(value).strip())
    return _stable_unique([part for part in parts if not is_missing(part)])


def normalize_inheritance_model(value: Any) -> str:
    model = str(value).strip()
    return INHERITANCE_ALIASES.get(model, model)


def _inheritance_pair_score(candidate: str, expected: str) -> float:
    if candidate == expected:
        if candidate == "AR_compound_heterozygous":
            return 80.0
        return 100.0

    return 0.0


def calculate_inheritance_score(
    inheritance_category: Any,
    expected_inheritance: Sequence[Any],
    *,
    ad_from_compound_pair: bool = False,
) -> float:
    """Score only categorical case-inheritance compatibility."""
    if not expected_inheritance:
        return 50.0

    candidates = [
        normalize_inheritance_model(model)
        for model in parse_inheritance_categories(inheritance_category)
    ]
    if not candidates:
        return 50.0

    expected = [
        normalize_inheritance_model(model)
        for model in expected_inheritance
        if not is_missing(model)
    ]
    if not expected:
        return 50.0

    direct_score = max(
        _inheritance_pair_score(candidate, expected_model)
        for candidate in candidates
        for expected_model in expected
    )
    if direct_score > 0:
        return direct_score

    expected_ar = any(
        model in {"AR_homozygous", "AR_compound_heterozygous"}
        for model in expected
    )
    selected_ad = any(model in {"AD", "AD_de_novo"} for model in candidates)
    if ad_from_compound_pair and expected_ar and selected_ad:
        return 25.0
    return 0.0


def effective_population_af(row: Mapping[str, Any]) -> float | None:
    """Return the Stage 3 filtering AF: max(global AF, grpmax AF)."""
    values = [
        value
        for value in (
            to_float(row.get("gnomAD_AF")),
            to_float(row.get("gnomAD_grpmax_AF")),
        )
        if value is not None
    ]
    if not values:
        return None
    return max(values)


def score_af_rarity(
    effective_af: float | None,
    threshold: float,
) -> float:
    """Score AF on the approved 0-10 logarithmic rarity scale."""
    if effective_af is None:
        return 5.0
    if threshold <= 0:
        raise ValueError("AF rarity threshold must be positive.")
    if effective_af <= 0:
        return 10.0

    log_ratio = math.log10(threshold / effective_af)
    return linear_score(
        log_ratio,
        RARITY_LOG_RATIO_ANCHORS,
        missing_score=5.0,
    )


def score_heterozygous_count(value: Any) -> float:
    count = to_float(value)
    if count is None:
        return 5.0
    if count < 0:
        raise Stage6InputError(f"Negative gnomAD heterozygote count: {value!r}")
    if count == 0:
        return 10.0
    if count <= 5:
        return 8.0
    if count <= 15:
        return 6.0
    if count <= 30:
        return 4.0
    if count < 50:
        return 2.0
    return 0.0


def nuclear_population_counts_are_zero(row: Mapping[str, Any]) -> bool:
    """Return True only when all three nuclear counts are present and zero."""
    fields = ("gnomAD_het", "gnomAD_hom", "gnomAD_hemi")
    if any(is_missing(row.get(field)) for field in fields):
        return False
    return all(to_float(row.get(field)) == 0.0 for field in fields)


def score_row_af_rarity(
    row: Mapping[str, Any],
    threshold: float,
    inheritance_model: str | None,
) -> float:
    """Score AF, treating explicit 0/0/0 nuclear counts as population absence."""
    effective_af = effective_population_af(row)
    if (
        effective_af is None
        and inheritance_model != "Mitochondrial"
        and nuclear_population_counts_are_zero(row)
    ):
        effective_af = 0.0
    return score_af_rarity(effective_af, threshold)


def score_homozygous_count(value: Any) -> float:
    count = to_float(value)
    if count is None:
        return 5.0
    if count < 0:
        raise Stage6InputError(f"Negative gnomAD homozygote count: {value!r}")
    if count == 0:
        return 10.0
    if count == 1:
        return 6.0
    if count == 2:
        return 2.0
    return 0.0


def score_compound_allele_count(row: Mapping[str, Any]) -> float:
    """Return the approved mean heterozygous/homozygous support per allele."""
    return (
        score_heterozygous_count(row.get("gnomAD_het"))
        + score_homozygous_count(row.get("gnomAD_hom"))
    ) / 2.0


def parse_compound_genes(value: Any) -> list[str]:
    if is_missing(value):
        return []
    return _stable_unique(
        [part for part in split_aligned(value) if not is_missing(part)]
    )


def parse_compound_partner_keys(
    value: Any,
    variant_index: Mapping[str, int],
) -> list[str]:
    """Parse Stage 3 partner keys while preserving the key's internal pipe."""
    if is_missing(value):
        return []

    raw_keys = re.split(r"\s+\|\s+", str(value).strip())
    partner_keys = _stable_unique(
        [key.strip() for key in raw_keys if not is_missing(key)]
    )
    unknown = [key for key in partner_keys if key not in variant_index]
    if unknown:
        raise Stage6InputError(
            "Compound_het_partner references variants absent from the Stage 5 "
            f"table: {', '.join(unknown)}"
        )
    return partner_keys


def validate_compound_partners_for_gene(
    row: Mapping[str, Any],
    gene: str,
    records: Sequence[Mapping[str, Any]],
    variant_index: Mapping[str, int],
) -> list[Mapping[str, Any]]:
    """Return reciprocal, different-locus partners for one compound gene."""
    row_key = variant_key(row)
    compound_genes = {
        item.casefold(): item for item in parse_compound_genes(row.get("Compound_het_gene"))
    }
    if gene.casefold() not in compound_genes:
        return []

    partner_keys = parse_compound_partner_keys(
        row.get("Compound_het_partner"),
        variant_index,
    )
    if not partner_keys:
        raise Stage6InputError(
            f"Compound-het candidate {row_key} has no listed partner for {gene}."
        )

    partners: list[Mapping[str, Any]] = []
    for partner_key in partner_keys:
        if partner_key == row_key:
            raise Stage6InputError(
                f"Compound-het candidate {row_key} lists itself as partner for {gene}."
            )
        partner = records[variant_index[partner_key]]
        partner_genes = {
            item.casefold()
            for item in parse_compound_genes(partner.get("Compound_het_gene"))
        }
        if gene.casefold() not in partner_genes:
            continue
        reciprocal_keys = parse_compound_partner_keys(
            partner.get("Compound_het_partner"),
            variant_index,
        )
        if row_key not in reciprocal_keys:
            raise Stage6InputError(
                f"Compound partner {partner_key} is not reciprocal with {row_key} "
                f"for {gene}."
            )
        if str(partner.get("Chr_Pos", "")).strip() == str(
            row.get("Chr_Pos", "")
        ).strip():
            raise Stage6InputError(
                f"Compound pair {row_key} / {partner_key} uses the same locus."
            )
        partners.append(partner)

    if not partners:
        raise Stage6InputError(
            f"Compound-het candidate {row_key} has no same-gene partner for {gene}."
        )

    return partners


def variant_key(row: Mapping[str, Any]) -> str:
    """Return the canonical in-memory key for one retained variant."""
    chromosome_position = str(row.get("Chr_Pos", "")).strip()
    ref_alt = str(row.get("Variant_REF_ALT", "")).strip()
    return f"{chromosome_position}|{ref_alt}"


def build_variant_index(frame: pd.DataFrame) -> dict[str, int]:
    """Map each retained variant key to its zero-based input row position."""
    index: dict[str, int] = {}
    key_columns = frame[["Chr_Pos", "Variant_REF_ALT"]]
    for position, (chromosome_position, ref_alt) in enumerate(
        key_columns.itertuples(index=False, name=None)
    ):
        key = f"{str(chromosome_position).strip()}|{str(ref_alt).strip()}"
        if key in index:
            raise Stage6InputError(f"Duplicate variant key in index: {key}")
        index[key] = position
    return index


def _normalized_text(value: str) -> str:
    """Normalize spacing only for safe repeated-value comparisons."""
    return re.sub(r"\s+", " ", value).strip()


def build_gene_cache(
    frame: pd.DataFrame,
    context: Mapping[str, Any],
) -> dict[str, dict[str, Any]]:
    """Build reusable, case-specific gene evidence without row expansion."""
    explicit_case_genes = {
        str(gene).strip().casefold()
        for field in ("case_genes", "manual_genes")
        for gene in context.get(field, [])
        if not is_missing(gene)
    }

    cache: dict[str, dict[str, Any]] = {}
    value_sources: dict[tuple[str, str], str] = {}

    for row in frame.to_dict(orient="records"):
        genes = parse_row_genes(row)
        if not genes:
            continue

        row_case_genes: set[str] = set()
        case_match = row.get("Case_gene_match")
        parsed_case_match = to_bool(case_match)
        if parsed_case_match is True and len(genes) == 1:
            row_case_genes.add(genes[0].casefold())
        elif parsed_case_match is None and not is_missing(case_match):
            case_match_text = str(case_match)
            row_case_genes.update(
                gene.casefold()
                for gene in genes
                if re.search(
                    rf"(?<![A-Za-z0-9_-]){re.escape(gene)}(?![A-Za-z0-9_-])",
                    case_match_text,
                    flags=re.I,
                )
            )

        parsed_fields = {
            field: parse_gene_blocks(row.get(field), genes, field_name=field)
            for field in GENE_CACHE_FIELDS
        }
        row_key = variant_key(row)

        for gene in genes:
            gene_key = gene.casefold()
            is_case_gene = (
                gene_key in explicit_case_genes or gene_key in row_case_genes
            )
            entry = cache.setdefault(
                gene,
                {
                    **{field: None for field in GENE_CACHE_FIELDS},
                    "is_case_gene": is_case_gene,
                },
            )
            entry["is_case_gene"] = bool(entry["is_case_gene"]) or is_case_gene

            for field in GENE_CACHE_FIELDS:
                parsed_value = parsed_fields[field][gene]
                if parsed_value is None:
                    continue

                if field in GENE_CACHE_NUMERIC_FIELDS:
                    new_value: Any = to_float(parsed_value)
                    if new_value is None:
                        raise Stage6InputError(
                            f"{field} is not numeric for {gene} at {row_key}: "
                            f"{parsed_value!r}"
                        )
                else:
                    new_value = parsed_value

                old_value = entry[field]
                if old_value is None:
                    entry[field] = new_value
                    value_sources[(gene, field)] = row_key
                    continue

                values_match = old_value == new_value
                if isinstance(old_value, str) and isinstance(new_value, str):
                    values_match = _normalized_text(old_value) == _normalized_text(
                        new_value
                    )

                if not values_match:
                    first_key = value_sources[(gene, field)]
                    raise Stage6InputError(
                        f"Inconsistent gene-level {field} for {gene}: "
                        f"{first_key} has {old_value!r}, while "
                        f"{row_key} has {new_value!r}."
                    )

    for entry in cache.values():
        entry["case_relevance_score"] = calculate_case_relevance_score(
            entry,
            context,
        )

    return cache


def _score_route(
    route: str,
    row: Mapping[str, Any],
    transcript: str,
    gene: str | None,
    transcript_evidence: Mapping[str, Mapping[str, str | None]],
    mane_transcripts: Sequence[str],
    gene_cache: Mapping[str, Mapping[str, Any]],
    warnings: list[str],
) -> tuple[float, float]:
    """Return route points out of 70 and observed route weight."""
    route = route.upper()

    if route == "COMMON":
        return 0.0, 0.0

    if route == "MISSENSE":
        evidence = transcript_evidence["Missense_evidence"].get(transcript)
        revel = extract_labeled_float(evidence, "REVEL")
        alpha_missense = extract_labeled_float(evidence, "AlphaMissense")
        protein_context = transcript_evidence["Protein_context"].get(transcript)
        protein_score, protein_observed = score_protein_context(protein_context)
        points = (
            linear_score(revel, REVEL_ANCHORS) * 3.0
            + linear_score(alpha_missense, ALPHAMISSENSE_ANCHORS) * 3.0
            + protein_score
        )
        observed_weight = (
            (30.0 if revel is not None else 0.0)
            + (30.0 if alpha_missense is not None else 0.0)
            + (10.0 if protein_observed else 0.0)
        )
        return points, observed_weight

    if route == "SPLICE":
        prediction = transcript_evidence["Splice_prediction"].get(transcript)
        splice_ai = extract_labeled_float(prediction, "SpliceAI")
        splice_context = transcript_evidence["Splice_context"].get(transcript)
        context_score, context_observed = score_splice_distance(splice_context)
        if not context_observed:
            warnings.append(
                f"{variant_key(row)}|{transcript}: SPLICE route has no "
                "transcript-specific splice distance"
            )
        points = linear_score(splice_ai, SPLICEAI_ANCHORS) * 4.0 + context_score * 3.0
        observed_weight = (40.0 if splice_ai is not None else 0.0) + (
            30.0 if context_observed else 0.0
        )
        return points, observed_weight

    if route == "LOF":
        lof_value = transcript_evidence["LoF_evidence"].get(transcript)
        nmd_value = transcript_evidence["NMD_context"].get(transcript)
        lof_score, lof_observed = score_lof_context(lof_value, nmd_value)
        loeuf = gene_cache.get(gene, {}).get("LOEUF") if gene else None
        loeuf_score, loeuf_observed = score_loeuf(loeuf)
        if mane_transcripts:
            mane_score = 10.0 if transcript in mane_transcripts else 5.0
            mane_observed = True
        else:
            mane_score = 5.0
            mane_observed = False
        points = lof_score * 4.0 + loeuf_score * 2.0 + mane_score
        observed_weight = (
            (40.0 if lof_observed else 0.0)
            + (20.0 if loeuf_observed else 0.0)
            + (10.0 if mane_observed else 0.0)
        )
        return points, observed_weight

    if route in {"INFRAME", "INFRAME_INDEL", "IN-FRAME"}:
        inframe_value = transcript_evidence["Inframe_indel_evidence"].get(
            transcript
        )
        functional_score, functional_observed = score_inframe_evidence(
            inframe_value
        )
        protein_value = transcript_evidence["Protein_context"].get(transcript)
        protein_score, protein_observed = score_protein_context(protein_value)
        points = functional_score * 5.0 + protein_score * 2.0
        observed_weight = (50.0 if functional_observed else 0.0) + (
            20.0 if protein_observed else 0.0
        )
        return points, observed_weight

    if route == "REGULATORY":
        annotation_score, annotation_observed = score_regulatory_annotation(
            row.get("Regulatory_context"),
            row.get("Regulatory_relevant"),
        )
        tss_value = transcript_evidence["TSS_context"].get(transcript)
        tss_score, tss_observed = score_tss_distance(tss_value)
        remm_score, remm_observed = score_remm(row.get("Regulatory_prediction"))
        points = annotation_score * 3.0 + tss_score * 2.0 + remm_score * 2.0
        observed_weight = (
            (30.0 if annotation_observed else 0.0)
            + (20.0 if tss_observed else 0.0)
            + (20.0 if remm_observed else 0.0)
        )
        return points, observed_weight

    raise Stage6InputError(
        f"Unsupported pathogenicity route {route!r} at {variant_key(row)}."
    )


def score_transcript_pathogenicity(
    row: Mapping[str, Any],
    gene_cache: Mapping[str, Mapping[str, Any]],
    warnings: list[str] | None = None,
) -> list[TranscriptPathogenicity]:
    """Calculate the winning pathogenicity route for every row transcript."""
    local_warnings = warnings if warnings is not None else []
    mitochondrial = is_mitochondrial_variant(row)
    if mitochondrial:
        local_warnings.append(
            f"{variant_key(row)}: partial mitochondrial scoring; only existing "
            "COMMON pathogenicity evidence is scored. MT-specific evidence and "
            "rarity are unavailable; rarity uses the existing missing default (50/100)."
        )

    transcript_gene_map = parse_transcript_gene_map(row)
    transcripts = list(transcript_gene_map)
    transcript_routes = parse_transcript_routes(
        row.get("Transcript_routes"),
        transcripts,
    )
    row_routes = parse_routes(row.get("Evidence_routes"))

    transcript_evidence = {
        field: parse_tx_blocks(
            row.get(field),
            transcripts=transcripts,
            warnings=local_warnings,
            field_name=field,
        )
        for field in TRANSCRIPT_EVIDENCE_FIELDS
    }
    mane_transcripts = parse_mane_transcripts(
        row.get("MANE_Transcript_ID"),
        transcripts,
    )

    cadd_points, cadd_observed = score_cadd(row.get("CADD"))
    conservation_points, conservation_sources = score_conservation(
        row.get("phyloP100way"),
        row.get("phastCons100way"),
    )
    common_points = cadd_points + conservation_points
    common_observed_weight = (20.0 if cadd_observed else 0.0) + (
        10.0 if conservation_sources else 0.0
    )
    clinvar_adjustment = clinvar_pathogenicity_adjustment(
        row.get("ClinVar_significance")
    )

    results: list[TranscriptPathogenicity] = []
    for transcript in transcripts:
        routes = [
            route
            for route in transcript_routes[transcript]
            if route != "COMMON"
        ]
        if not routes and len(transcripts) == 1:
            fallback_routes = [
                route
                for route in row_routes
                if route not in {"COMMON", "MITOCHONDRIAL"}
            ]
            if fallback_routes:
                routes.extend(fallback_routes)
                local_warnings.append(
                    f"{variant_key(row)}|{transcript}: using Evidence_routes "
                    "because Transcript_routes is missing"
                )
        if "REGULATORY" in row_routes and "REGULATORY" not in routes:
            routes.append("REGULATORY")
        if mitochondrial:
            routes = ["COMMON"]
        if not routes:
            routes = ["COMMON"]

        best_result: TranscriptPathogenicity | None = None
        for route in routes:
            route_points, route_observed_weight = _score_route(
                route=route,
                row=row,
                transcript=transcript,
                gene=transcript_gene_map[transcript],
                transcript_evidence=transcript_evidence,
                mane_transcripts=mane_transcripts,
                gene_cache=gene_cache,
                warnings=local_warnings,
            )
            result = TranscriptPathogenicity(
                transcript=transcript,
                gene=transcript_gene_map[transcript],
                route=route,
                score=clamp_score(
                    common_points + route_points + clinvar_adjustment
                ),
                coverage=(common_observed_weight + route_observed_weight) / 100.0,
            )
            if best_result is None or result.score > best_result.score:
                best_result = result

        assert best_result is not None
        results.append(best_result)

    return results


def select_best_transcript_hypotheses(
    results: Sequence[TranscriptPathogenicity],
) -> dict[str | None, TranscriptPathogenicity]:
    """Keep one transcript per gene plus one transcript-only hypothesis."""
    best_by_gene: dict[str | None, TranscriptPathogenicity] = {}
    for result in results:
        current = best_by_gene.get(result.gene)
        if current is None or result.score > current.score:
            best_by_gene[result.gene] = result
    return best_by_gene


def _normalized_match_text(value: Any) -> str:
    if is_missing(value):
        return ""
    return re.sub(r"[^a-z0-9]+", " ", str(value).casefold()).strip()


def _contains_case_term(evidence: Any, terms: Sequence[Any]) -> bool:
    normalized_evidence = _normalized_match_text(evidence)
    if not normalized_evidence:
        return False
    return any(
        normalized_term in normalized_evidence
        for term in terms
        if (normalized_term := _normalized_match_text(term))
    )


def extract_hpo_ids(value: Any) -> set[str]:
    if is_missing(value):
        return set()
    return {
        match.upper().replace("_", ":")
        for match in re.findall(r"HP[:_]\d{7}", str(value), flags=re.I)
    }


def calculate_case_relevance_score(
    gene_evidence: Mapping[str, Any],
    context: Mapping[str, Any],
) -> float:
    """Apply the approved hierarchical gene-to-case relevance rules."""
    if bool(gene_evidence.get("is_case_gene")):
        return 100.0

    disease_evidence = gene_evidence.get("Disease_associations")
    if _contains_case_term(disease_evidence, context.get("disease_names", [])):
        return 100.0

    applicable_scores: list[float] = []
    broad_domains = context.get("broad_disease_domains", [])
    if broad_domains and not is_missing(disease_evidence):
        if _contains_case_term(disease_evidence, broad_domains):
            applicable_scores.append(70.0)

    patient_hpo = extract_hpo_ids(" | ".join(map(str, context.get("hpo_terms", []))))
    gene_hpo = extract_hpo_ids(gene_evidence.get("Gene_HPO_terms"))
    if patient_hpo and gene_hpo:
        matched_hpo = extract_hpo_ids(gene_evidence.get("Matched_HPO_terms"))
        overlap = matched_hpo or (patient_hpo & gene_hpo)
        applicable_scores.append(100.0 * len(overlap & gene_hpo) / len(gene_hpo))

    if not applicable_scores:
        has_case_context = any(
            context.get(field, [])
            for field in (
                "disease_names",
                "hpo_terms",
                "broad_disease_domains",
                "manual_genes",
                "case_genes",
            )
        )
        return 20.0 if has_case_context else 50.0
    return clamp_score(max(applicable_scores))


def score_gencc_strength(value: Any) -> float:
    """Return the GenCC contribution out of 50 points."""
    if is_missing(value):
        return 25.0
    text = str(value).casefold()
    if re.search(r"\b(definitive|strong)\b", text):
        return 50.0
    if re.search(r"\b(moderate|supportive|limited|animal model only)\b", text):
        return 25.0
    if re.search(
        r"\b(disputed|refuted|no known disease relationship)\b",
        text,
    ):
        return 0.0
    return 25.0


def _inheritance_family(model: str | None) -> str | None:
    if model in {"AD", "AD_de_novo", "Mosaic"}:
        return "AD"
    if model in {"AR_homozygous", "AR_compound_heterozygous"}:
        return "AR"
    if model == "X_linked":
        return "X"
    if model == "Mitochondrial":
        return "MT"
    return None


def _evidence_inheritance_families(value: Any) -> set[str]:
    if is_missing(value):
        return set()
    text = str(value).casefold()
    families: set[str] = set()
    if "autosomal dominant" in text or re.search(r"\bAD\b", str(value)):
        families.add("AD")
    if "autosomal recessive" in text or re.search(r"\bAR\b", str(value)):
        families.add("AR")
    if "x-linked" in text or "x linked" in text:
        families.add("X")
    if "mitochond" in text or "maternal" in text:
        families.add("MT")
    return families


def score_moi_compatibility(model: str | None, value: Any) -> float:
    """Return the gene-disease MOI contribution out of 15 points."""
    candidate_family = _inheritance_family(model)
    evidence_families = _evidence_inheritance_families(value)
    if candidate_family is None or not evidence_families:
        return 7.5
    if model == "Mosaic" and "AD" in evidence_families:
        return 7.5
    if candidate_family in evidence_families:
        return 15.0
    return 0.0


def _route_landscape_terms(route: str) -> tuple[str, ...]:
    normalized = route.upper()
    if normalized == "MISSENSE":
        return ("missense",)
    if normalized == "SPLICE":
        return ("splice",)
    if normalized == "LOF":
        return ("frameshift", "nonsense", "start-loss", "start loss", "stop-loss", "stop loss")
    if normalized in {"INFRAME", "INFRAME_INDEL", "IN-FRAME"}:
        return ("in-frame", "inframe")
    if normalized == "REGULATORY":
        return ("regulatory", "intronic", "utr", "noncoding", "non-coding", "near-gene")
    return ()


def score_clinvar_class_match(
    route: str,
    landscape: Any,
    pathogenic_count: float | None,
) -> float:
    """Return the route/class match contribution out of 10 points."""
    terms = _route_landscape_terms(route)
    if not terms:
        return 5.0
    if is_missing(landscape):
        return 5.0 if pathogenic_count is None or pathogenic_count > 0 else 0.0

    text = str(landscape).casefold()
    for term in terms:
        match = re.search(rf"{re.escape(term)}\s*=\s*({NUMBER_PATTERN})", text)
        if match and (to_float(match.group(1)) or 0.0) > 0:
            return 10.0
    if pathogenic_count is None:
        return 5.0
    return 5.0 if pathogenic_count > 0 else 0.0


def score_clinvar_pathogenic_count(value: float | None) -> float:
    if value is None:
        return 2.5
    if value <= 0:
        return 0.0
    if value == 1:
        return 1.0
    if value <= 4:
        return 2.0
    if value <= 10:
        return 3.5
    return 5.0


def score_clinvar_pathogenic_fraction(value: float | None) -> float:
    if value is None:
        return 2.5
    if value <= 0:
        return 0.0
    if value < 0.05:
        return 1.0
    if value < 0.10:
        return 2.5
    return 5.0


def score_literature_support(value: Any) -> float:
    """Return the row-specific literature contribution out of 15 points."""
    if is_missing(value):
        return 7.5
    text = str(value).casefold()
    human_match = re.search(r"human-cases\s*=\s*(\d+)", text)
    human_cases = int(human_match.group(1)) if human_match else 0
    if human_cases >= 2:
        return 15.0
    if human_cases == 1:
        return 12.0
    functional_counts = re.findall(r"(?:functional|animal)\s*=\s*(\d+)", text)
    if any(int(count) > 0 for count in functional_counts):
        return 7.5
    reviewed_match = re.search(r"total publications reviewed\s*:\s*(\d+)", text)
    if (reviewed_match and int(reviewed_match.group(1)) == 0) or "no publications" in text:
        return 0.0
    return 7.5


def calculate_gene_disease_score(
    gene: str | None,
    pathogenicity_route: str,
    inheritance_model: str | None,
    gene_cache: Mapping[str, Mapping[str, Any]],
    literature_support: Any,
) -> float:
    if gene is None:
        return NO_GENE_GENE_DISEASE_SCORE
    evidence = gene_cache.get(gene, {})
    count = evidence.get("ClinVar_pathogenic_count")
    fraction = evidence.get("ClinVar_pathogenic_fraction")
    score = (
        score_gencc_strength(evidence.get("Disease_associations"))
        + score_moi_compatibility(
            inheritance_model,
            evidence.get("Inheritance_and_mechanism"),
        )
        + score_clinvar_class_match(
            pathogenicity_route,
            evidence.get("Known_variant_landscape"),
            count,
        )
        + score_clinvar_pathogenic_count(count)
        + score_clinvar_pathogenic_fraction(fraction)
        + score_literature_support(literature_support)
    )
    return clamp_score(score)


def select_singleton_inheritance_model(
    inheritance_category: Any,
    expected_inheritance: Sequence[Any],
) -> str | None:
    """Choose one non-compound route without constructing route combinations."""
    candidates = []
    for model in parse_inheritance_categories(inheritance_category):
        normalized_model = normalize_inheritance_model(model)
        if normalized_model != "AR_compound_heterozygous":
            candidates.append(normalized_model)
    if not candidates:
        return None
    if not expected_inheritance:
        return candidates[0]
    expected = [
        normalize_inheritance_model(model)
        for model in expected_inheritance
        if not is_missing(model)
    ]
    if not expected:
        return candidates[0]
    return max(
        enumerate(candidates),
        key=lambda item: (
            max(_inheritance_pair_score(item[1], model) for model in expected),
            -item[0],
        ),
    )[1]


def calculate_singleton_rarity_score(
    row: Mapping[str, Any],
    inheritance_model: str | None,
    is_case_gene: bool,
) -> float:
    model = inheritance_model or "AD"
    if is_mitochondrial_variant(row) or model == "Mitochondrial":
        # MT AF/count definitions do not support the nuclear rarity model.
        af_score = score_af_rarity(None, RARITY_AF_THRESHOLDS["Mitochondrial"])
        count_score = score_homozygous_count(None)
        return clamp_score((0.8 * af_score + 0.2 * count_score) * 10.0)
    if model in {"AD", "AD_de_novo"}:
        threshold_key = "AD_case_gene" if is_case_gene else "AD"
        count_score = score_heterozygous_count(row.get("gnomAD_het"))
    elif model == "AR_homozygous":
        threshold_key = "AR_homozygous"
        count_score = score_homozygous_count(row.get("gnomAD_hom"))
    elif model == "X_linked":
        threshold_key = "X_linked"
        count_score = score_homozygous_count(row.get("gnomAD_hemi"))
    elif model == "Mosaic":
        threshold_key = "Mosaic"
        count_score = score_heterozygous_count(row.get("gnomAD_het"))
    else:
        threshold_key = "AD"
        count_score = score_heterozygous_count(row.get("gnomAD_het"))
    af_score = score_row_af_rarity(
        row,
        RARITY_AF_THRESHOLDS[threshold_key],
        model,
    )
    return clamp_score((0.8 * af_score + 0.2 * count_score) * 10.0)


def calculate_compound_allele_rarity_score(
    row: Mapping[str, Any],
    is_case_gene: bool,
) -> float:
    threshold_key = (
        "AR_compound_heterozygous_case_gene"
        if is_case_gene
        else "AR_compound_heterozygous"
    )
    threshold = RARITY_AF_THRESHOLDS[threshold_key]
    af_support = score_row_af_rarity(
        row,
        threshold,
        "AR_compound_heterozygous",
    )
    count_support = score_compound_allele_count(row)
    return clamp_score((0.8 * af_support + 0.2 * count_support) * 10.0)


def calculate_final_score(
    pathogenicity: float,
    inheritance: float,
    gene_disease: float,
    case_relevance: float,
    rarity: float,
) -> float:
    return round(
        0.30 * pathogenicity
        + 0.20 * inheritance
        + 0.25 * gene_disease
        + 0.15 * case_relevance
        + 0.10 * rarity,
        2,
    )


def _candidate_choice_key(candidate: Candidate) -> tuple[Any, ...]:
    return (
        candidate.has_priority_exception,
        candidate.final_score,
        candidate.variant_pathogenicity_score,
        candidate.gene_disease_evidence_score,
        -candidate.gene_order,
    )


def choose_best_candidate(candidates: Sequence[Candidate]) -> Candidate:
    if not candidates:
        raise Stage6InputError("A variant produced no scoreable candidate.")
    return max(candidates, key=_candidate_choice_key)


def build_singleton_candidates(
    row: Mapping[str, Any],
    gene_cache: Mapping[str, Mapping[str, Any]],
    context: Mapping[str, Any],
    warnings: list[str],
) -> list[Candidate]:
    pathogenicity = select_best_transcript_hypotheses(
        score_transcript_pathogenicity(row, gene_cache, warnings)
    )
    row_genes = parse_row_genes(row)
    literature_by_gene = parse_gene_blocks(
        row.get("Literature_support"),
        row_genes,
        field_name="Literature_support",
    )
    gene_order = {gene: order for order, gene in enumerate(row_genes)}
    gene_order[None] = len(row_genes)
    expected_inheritance = context.get("expected_inheritance", [])
    inheritance_model = select_singleton_inheritance_model(
        row.get("inheritance_category"),
        expected_inheritance,
    )
    inheritance_score = calculate_inheritance_score(
        inheritance_model or "",
        expected_inheritance,
        ad_from_compound_pair=(
            inheritance_model in {"AD", "AD_de_novo"}
            and bool(parse_compound_genes(row.get("Compound_het_gene")))
        ),
    )
    is_pathogenic = clinvar_has_pathogenic_evidence(
        row.get("ClinVar_significance")
    )

    candidates: list[Candidate] = []
    for gene, transcript_result in pathogenicity.items():
        is_case_gene = bool(gene and gene_cache.get(gene, {}).get("is_case_gene"))
        gene_disease = round(
            calculate_gene_disease_score(
                gene,
                transcript_result.route,
                inheritance_model,
                gene_cache,
                literature_by_gene.get(gene),
            ),
            2,
        )
        case_relevance = round(
            float(gene_cache.get(gene, {}).get("case_relevance_score", 0.0))
            if gene
            else NO_GENE_CASE_RELEVANCE_SCORE,
            2,
        )
        rarity = round(
            calculate_singleton_rarity_score(row, inheritance_model, is_case_gene),
            2,
        )
        pathogenicity_score = round(transcript_result.score, 2)
        candidates.append(
            Candidate(
                gene=gene,
                gene_order=gene_order.get(gene, len(row_genes)),
                transcript=transcript_result.transcript,
                pathogenicity_route=transcript_result.route,
                variant_pathogenicity_score=pathogenicity_score,
                pathogenicity_evidence_coverage=round(
                    transcript_result.coverage,
                    4,
                ),
                inheritance_score=round(inheritance_score, 2),
                gene_disease_evidence_score=gene_disease,
                case_relevance_score=case_relevance,
                rarity_score=rarity,
                final_score=calculate_final_score(
                    pathogenicity_score,
                    inheritance_score,
                    gene_disease,
                    case_relevance,
                    rarity,
                ),
                is_case_gene=is_case_gene,
                is_pathogenic_clinvar=is_pathogenic,
                inheritance_model=inheritance_model or "",
            )
        )
    if "AR_compound_heterozygous" in parse_inheritance_categories(row.get("inheritance_category")):
        compound_genes = {gene.casefold() for gene in parse_compound_genes(row.get("Compound_het_gene"))}
        # Compare fully scored hypotheses before choosing the gene. An additional
        # allele in an AR gene stays AR even when another pair represents that
        # gene in the report. AD is the alternative hypothesis for other genes.
        ad_candidates = build_singleton_candidates(
            dict(row, inheritance_category="AD"), gene_cache, context, [],
        )
        ad_by_gene = {candidate.gene: candidate for candidate in ad_candidates}
        candidates = [
            _make_compound_candidate(candidate, row, "", gene_cache, context)
            if candidate.gene and candidate.gene.casefold() in compound_genes
            else ad_by_gene[candidate.gene] if candidate.gene else candidate
            for candidate in candidates
        ]
    return candidates


def _candidate_for_gene(
    candidates: Sequence[Candidate],
    gene: str,
) -> Candidate:
    for candidate in candidates:
        if candidate.gene and candidate.gene.casefold() == gene.casefold():
            return candidate
    raise Stage6InputError(
        f"Compound gene {gene!r} has no transcript-aligned scoring hypothesis."
    )


def _make_compound_candidate(
    base: Candidate,
    row: Mapping[str, Any],
    partner_key: str,
    gene_cache: Mapping[str, Mapping[str, Any]],
    context: Mapping[str, Any],
) -> Candidate:
    inheritance_score = calculate_inheritance_score(
        "AR_compound_heterozygous",
        context.get("expected_inheritance", []),
    )
    gene_disease = round(
        calculate_gene_disease_score(
            base.gene,
            base.pathogenicity_route,
            "AR_compound_heterozygous",
            gene_cache,
            parse_gene_blocks(
                row.get("Literature_support"),
                parse_row_genes(row),
                field_name="Literature_support",
            ).get(base.gene),
        ),
        2,
    )
    is_case_gene = bool(
        base.gene and gene_cache.get(base.gene, {}).get("is_case_gene")
    )
    allele_rarity = calculate_compound_allele_rarity_score(row, is_case_gene)
    final_score = calculate_final_score(
        base.variant_pathogenicity_score,
        inheritance_score,
        gene_disease,
        base.case_relevance_score,
        allele_rarity,
    )
    return base._replace(
        inheritance_score=round(inheritance_score, 2),
        gene_disease_evidence_score=gene_disease,
        rarity_score=round(allele_rarity, 2),
        final_score=final_score,
        scored_compound_partner=partner_key,
        inheritance_model="AR_compound_heterozygous",
    )


def build_compound_pairs(
    records: Sequence[Mapping[str, Any]],
    singleton_candidates: Sequence[Sequence[Candidate]],
    variant_index: Mapping[str, int],
    gene_cache: Mapping[str, Mapping[str, Any]],
    context: Mapping[str, Any],
) -> list[CompoundPair]:
    pairs: list[CompoundPair] = []
    seen: set[tuple[int, int, str]] = set()
    for left_index, row in enumerate(records):
        left_key = variant_key(row)
        for gene in parse_compound_genes(row.get("Compound_het_gene")):
            partners = validate_compound_partners_for_gene(
                row,
                gene,
                records,
                variant_index,
            )
            for partner in partners:
                right_key = variant_key(partner)
                right_index = variant_index[right_key]
                pair_key = (
                    min(left_index, right_index),
                    max(left_index, right_index),
                    gene.casefold(),
                )
                if pair_key in seen:
                    continue
                seen.add(pair_key)
                left_base = _candidate_for_gene(singleton_candidates[left_index], gene)
                right_base = _candidate_for_gene(singleton_candidates[right_index], gene)
                left_candidate = _make_compound_candidate(
                    left_base,
                    row,
                    right_key,
                    gene_cache,
                    context,
                )
                right_candidate = _make_compound_candidate(
                    right_base,
                    partner,
                    left_key,
                    gene_cache,
                    context,
                )
                pair_score = round(
                    (left_candidate.final_score + right_candidate.final_score) / 2.0,
                    2,
                )
                pairs.append(
                    CompoundPair(
                        gene=left_base.gene or gene,
                        left_index=left_index,
                        right_index=right_index,
                        left_candidate=left_candidate._replace(
                            compound_pair_score=pair_score,
                        ),
                        right_candidate=right_candidate._replace(
                            compound_pair_score=pair_score,
                        ),
                        pair_score=pair_score,
                    )
                )
    return pairs


def _compound_beats_singleton(
    compound_candidate: Candidate,
    pair_score: float,
    singleton_candidate: Candidate,
) -> bool:
    if (
        compound_candidate.gene and singleton_candidate.gene
        and compound_candidate.gene.casefold() == singleton_candidate.gene.casefold()
        and singleton_candidate.inheritance_model == "AR_compound_heterozygous"
    ):
        # Same-gene AR alleles support the pair rather than competing with its
        # mean score. Only an alternative gene hypothesis can displace it.
        return True
    if (
        compound_candidate.has_priority_exception
        != singleton_candidate.has_priority_exception
    ):
        return compound_candidate.has_priority_exception
    return pair_score >= singleton_candidate.final_score


def select_final_candidates(
    singleton_candidates: Sequence[Sequence[Candidate]],
    compound_pairs: Sequence[CompoundPair],
) -> list[Candidate]:
    selected = [choose_best_candidate(candidates) for candidates in singleton_candidates]
    occupied: set[int] = set()
    selected_genes: set[str] = set()
    ordered_pairs = sorted(
        compound_pairs,
        key=lambda pair: (
            pair.left_candidate.has_priority_exception
            + pair.right_candidate.has_priority_exception,
            pair.pair_score,
            pair.left_candidate.variant_pathogenicity_score
            + pair.right_candidate.variant_pathogenicity_score,
            -min(pair.left_index, pair.right_index),
            -max(pair.left_index, pair.right_index),
            pair.gene,
        ),
        reverse=True,
    )
    for pair in ordered_pairs:
        gene_key = pair.gene.casefold()
        if gene_key in selected_genes:
            continue
        if pair.left_index in occupied or pair.right_index in occupied:
            continue
        if not _compound_beats_singleton(
            pair.left_candidate,
            pair.pair_score,
            selected[pair.left_index],
        ):
            continue
        if not _compound_beats_singleton(
            pair.right_candidate,
            pair.pair_score,
            selected[pair.right_index],
        ):
            continue
        selected[pair.left_index] = pair.left_candidate
        selected[pair.right_index] = pair.right_candidate
        occupied.update((pair.left_index, pair.right_index))
        selected_genes.add(gene_key)

    # Additional alleles in a gene with a selected pair remain compound
    # candidates for manual review. They keep compound inheritance, MOI and
    # allele-specific rarity, but only the selected pair receives pair fields.
    for pair in ordered_pairs:
        if pair.gene.casefold() not in selected_genes:
            continue
        for row_index, compound_candidate in (
            (pair.left_index, pair.left_candidate),
            (pair.right_index, pair.right_candidate),
        ):
            if row_index in occupied:
                continue
            unpaired_candidate = compound_candidate._replace(
                scored_compound_partner="",
                compound_pair_score=None,
            )
            if _candidate_choice_key(unpaired_candidate) > _candidate_choice_key(
                selected[row_index]
            ):
                selected[row_index] = unpaired_candidate
    return selected


def recalculate_ad_fallbacks(records, selected, variant_index, gene_cache, context, warnings):
    """Mark AD only when the winning gene is outside the variant's AR genes."""
    result = list(selected)
    for index, (row, candidate) in enumerate(zip(records, selected, strict=True)):
        if "AR_compound_heterozygous" not in parse_inheritance_categories(row.get("inheritance_category")):
            continue
        compound_genes = {gene.casefold() for gene in parse_compound_genes(row.get("Compound_het_gene"))}
        if not candidate.gene or candidate.gene.casefold() in compound_genes:
            continue
        ad_row = dict(row, inheritance_category="AD")
        candidates = build_singleton_candidates(ad_row, gene_cache, context, warnings)
        ad = next((item for item in candidates if item.gene == candidate.gene), None)
        if ad is None:
            raise Stage6InputError("AD fallback lost the selected gene hypothesis")
        result[index] = ad._replace(ad_fallback=True, scored_compound_partner="", compound_pair_score=None)
    return result


def _candidate_output(candidate: Candidate) -> dict[str, Any]:
    return {
        "Variant_pathogenicity_score": candidate.variant_pathogenicity_score,
        "Pathogenicity_evidence_coverage": candidate.pathogenicity_evidence_coverage,
        "Inheritance_score": candidate.inheritance_score,
        "Gene_disease_evidence_score": candidate.gene_disease_evidence_score,
        "Case_relevance_score": candidate.case_relevance_score,
        "Rarity_score": candidate.rarity_score,
        "Scored_gene": candidate.gene or NO_GENE_SCORED_GENE,
        "Scored_transcript": candidate.transcript,
        "Scored_route": candidate.pathogenicity_route,
        "Scored_inheritance": candidate.inheritance_model,
        "AD_fallback": candidate.ad_fallback,
        "Scored_compound_partner": candidate.scored_compound_partner,
        "Compound_pair_score": (
            "" if candidate.compound_pair_score is None else candidate.compound_pair_score
        ),
        "Final_score": candidate.final_score,
    }


def rank_scored_output(
    output: pd.DataFrame,
    selected: Sequence[Candidate],
) -> pd.DataFrame:
    """Keep each scored gene together, with its selected pair first."""
    group_members: dict[tuple[str, Any], list[int]] = {}
    for index, candidate in enumerate(selected):
        group_key = (
            ("gene", candidate.gene.casefold())
            if candidate.gene
            else ("row", index)
        )
        group_members.setdefault(group_key, []).append(index)

    def group_sort_key(item: tuple[tuple[str, Any], list[int]]) -> tuple[Any, ...]:
        _, indices = item
        candidates = [selected[index] for index in indices]
        pair_scores = [
            candidate.compound_pair_score
            for candidate in candidates
            if candidate.compound_pair_score is not None
        ]
        ranking_value = (
            max(pair_scores)
            if pair_scores
            else max(candidate.final_score for candidate in candidates)
        )
        return (
            any(candidate.is_case_gene for candidate in candidates),
            ranking_value,
            max(candidate.final_score for candidate in candidates),
            max(candidate.variant_pathogenicity_score for candidate in candidates),
            max(candidate.gene_disease_evidence_score for candidate in candidates),
            -min(indices),
        )

    ordered_groups = sorted(
        group_members.items(),
        key=group_sort_key,
        reverse=True,
    )
    group_order = {
        group_key: order
        for order, (group_key, _) in enumerate(ordered_groups)
    }

    ranked = output.copy()
    ranked["__group_order"] = [
        group_order[
            ("gene", candidate.gene.casefold())
            if candidate.gene
            else ("row", index)
        ]
        for index, candidate in enumerate(selected)
    ]
    ranked["__selected_pair"] = [
        candidate.compound_pair_score is not None
        for candidate in selected
    ]
    ranked["__input_order"] = range(len(ranked))
    return (
        ranked.sort_values(
            by=[
                "__group_order",
                "__selected_pair",
                "Final_score",
                "Variant_pathogenicity_score",
                "Gene_disease_evidence_score",
                "__input_order",
            ],
            ascending=[True, False, False, False, False, True],
            kind="mergesort",
        )
        .drop(columns=["__group_order", "__selected_pair", "__input_order"])
        .reset_index(drop=True)
    )


def score_stage6(
    frame: pd.DataFrame,
    context: Mapping[str, Any],
    gene_cache: Mapping[str, Mapping[str, Any]],
    variant_index: Mapping[str, int],
    warnings: list[str],
    progress_every: int = 100,
) -> tuple[pd.DataFrame, list[Candidate]]:
    records = frame.to_dict(orient="records")
    singleton_candidates: list[list[Candidate]] = []
    total = len(records)
    for row_index, row in enumerate(records):
        singleton_candidates.append(
            build_singleton_candidates(
                row,
                gene_cache,
                context,
                warnings,
            )
        )
        if progress_every > 0 and (
            (row_index + 1) % progress_every == 0 or row_index + 1 == total
        ):
            print(f"Scoring variants: {row_index + 1:,} / {total:,}")

    compound_pairs = build_compound_pairs(
        records,
        singleton_candidates,
        variant_index,
        gene_cache,
        context,
    )
    selected = select_final_candidates(singleton_candidates, compound_pairs)
    selected = recalculate_ad_fallbacks(records, selected, variant_index, gene_cache, context, warnings)

    score_frame = pd.DataFrame([_candidate_output(candidate) for candidate in selected])
    output = pd.concat(
        [frame.reset_index(drop=True), score_frame],
        axis=1,
    )

    return rank_scored_output(output, selected), selected


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--case-context", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--progress-every", type=int, default=100)
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Load and validate inputs without producing a Stage 6 table.",
    )
    return parser.parse_args()


def args_from_snakemake(smk: Any) -> argparse.Namespace:
    """Translate the Snakemake object into the same arguments used by the CLI."""
    return argparse.Namespace(
        input=Path(str(smk.input.variants)),
        case_context=Path(str(smk.input.case_context)),
        output=Path(str(smk.output.variants)),
        progress_every=int(smk.params.progress_every),
        validate_only=False,
    )


def load_stage5(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise Stage6InputError(f"Stage 5 table not found: {path}")

    print("Loading Stage 5 input...")
    unicode_error: UnicodeDecodeError | None = None
    for encoding in ("utf-8-sig", "cp1252"):
        try:
            frame = pd.read_csv(
                path,
                dtype=str,
                keep_default_na=False,
                encoding=encoding,
            )
            if encoding != "utf-8-sig":
                print(
                    "Warning: Stage 5 input is not UTF-8; "
                    f"loaded using {encoding}."
                )
            return frame
        except UnicodeDecodeError as exc:
            unicode_error = exc
        except Exception as exc:
            raise Stage6InputError(f"Could not read Stage 5 table: {path}") from exc

    raise Stage6InputError(
        f"Could not decode Stage 5 table as UTF-8 or Windows-1252: {path}"
    ) from unicode_error


def load_case_context(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise Stage6InputError(f"Case Context not found: {path}")

    print("Loading Case Context...")
    try:
        with path.open("r", encoding="utf-8") as handle:
            context = yaml.safe_load(handle)
    except Exception as exc:
        raise Stage6InputError(f"Could not read Case Context: {path}") from exc

    if not isinstance(context, Mapping):
        raise Stage6InputError("Case Context must be a YAML mapping.")

    return dict(context)


def validate_stage5_structure(frame: pd.DataFrame) -> list[str]:
    warnings: list[str] = []
    validate_table(frame, "Stage 6 input", aligned=True)
    if frame.empty:
        raise Stage6InputError("Stage 5 table contains no variant rows.")
    required = REQUIRED_STAGE5_COLUMNS + REQUIRED_SCORING_COLUMNS
    missing = [name for name in required if name not in frame.columns]
    if missing:
        joined = ", ".join(missing)
        raise Stage6InputError(f"Stage 5 table is missing required columns: {joined}")

    duplicate_columns = frame.columns[frame.columns.duplicated()].tolist()
    if duplicate_columns:
        joined = ", ".join(duplicate_columns)
        raise Stage6InputError(f"Stage 5 table has duplicate columns: {joined}")

    duplicate_keys = frame.duplicated(subset=["Chr_Pos", "Variant_REF_ALT"], keep=False)
    if duplicate_keys.any():
        examples = (
            frame.loc[duplicate_keys, ["Chr_Pos", "Variant_REF_ALT"]]
            .head(5)
            .astype(str)
            .agg("|".join, axis=1)
            .tolist()
        )
        raise Stage6InputError(
            "Stage 5 table has duplicate variant keys. Examples: "
            + ", ".join(examples)
        )

    existing_output_columns = [
        name for name in STAGE6_OUTPUT_COLUMNS if name in frame.columns
    ]
    if existing_output_columns:
        joined = ", ".join(existing_output_columns)
        raise Stage6InputError(
            f"Stage 5 input already contains Stage 6 output columns: {joined}"
        )

    for row in frame.to_dict(orient="records"):
        transcript_gene_map = parse_transcript_gene_map(row)
        transcripts = list(transcript_gene_map)
        genes = parse_row_genes(row)

        if not genes:
            key = f"{row.get('Chr_Pos', '?')}|{row.get('Variant_REF_ALT', '?')}"
            warnings.append(
                f"{key}: no gene hypothesis; Scored_gene will be empty and "
                "gene-dependent scores will be 0"
            )

        for field_name in TRANSCRIPT_EVIDENCE_FIELDS:
            parse_tx_blocks(
                row.get(field_name),
                transcripts=transcripts,
                warnings=warnings,
                field_name=field_name,
            )

        for field_name in GENE_EVIDENCE_FIELDS:
            parse_gene_blocks(
                row.get(field_name),
                genes,
                warnings=warnings,
                field_name=field_name,
            )

        parse_transcript_routes(row.get("Transcript_routes"), transcripts)
        parse_mane_transcripts(row.get("MANE_Transcript_ID"), transcripts)

        for field_name in ("Nearest_splice_transcript", "Nearest_TSS_transcript"):
            raw_transcript = row.get(field_name)
            if is_missing(raw_transcript):
                continue
            if resolve_transcript_id(str(raw_transcript), transcripts) is None:
                key = f"{row.get('Chr_Pos', '?')}|{row.get('Variant_REF_ALT', '?')}"
                raise Stage6InputError(
                    f"{field_name} references an unknown transcript for {key}: "
                    f"{raw_transcript!r}."
                )

    return warnings


def validate_case_context(context: Mapping[str, Any]) -> None:
    missing = [name for name in REQUIRED_CASE_FIELDS if name not in context]
    if missing:
        joined = ", ".join(missing)
        raise Stage6InputError(f"Case Context is missing required fields: {joined}")

    malformed = [
        name
        for name in LIST_CASE_FIELDS
        if not isinstance(context[name], list)
    ]
    if malformed:
        joined = ", ".join(malformed)
        raise Stage6InputError(f"Case Context fields must be lists: {joined}")

    sex = str(context["sex"]).strip()
    if not sex:
        raise Stage6InputError("Case Context field 'sex' cannot be empty.")

    unknown_inheritance = [
        model
        for model in context["expected_inheritance"]
        if model not in EXPECTED_INHERITANCE_MODELS
    ]
    if unknown_inheritance:
        joined = ", ".join(map(str, unknown_inheritance))
        raise Stage6InputError(
            f"Case Context has unsupported expected inheritance models: {joined}"
        )


def validate_inputs(
    frame: pd.DataFrame,
    context: Mapping[str, Any],
) -> list[str]:
    print("Validating...")
    warnings = validate_stage5_structure(frame)
    validate_case_context(context)
    return warnings


def atomic_write_csv(frame: pd.DataFrame, destination: Path) -> None:
    """Write a CSV beside its destination and replace the target atomically."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None

    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            suffix=".tmp",
            prefix=f".{destination.name}.",
            dir=destination.parent,
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            frame.to_csv(handle, index=False)

        os.replace(temporary_path, destination)
    except Exception:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
        raise


def describe_scores(frame: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    """Return the required distribution summary for completed score columns."""
    rows: dict[str, dict[str, float]] = {}
    for column in columns:
        values = pd.to_numeric(frame[column], errors="coerce").dropna()
        rows[column] = {
            "count": float(values.count()),
            "mean": float(values.mean()),
            "median": float(values.median()),
            "min": float(values.min()),
            "p25": float(values.quantile(0.25)),
            "p75": float(values.quantile(0.75)),
            "max": float(values.max()),
        }
    return pd.DataFrame.from_dict(rows, orient="index")


def validate_scored_output(
    input_frame: pd.DataFrame,
    output_frame: pd.DataFrame,
) -> None:
    print("Validating final output...")
    validate_table(output_frame, "Stage 6 scored output", aligned=True)
    if len(input_frame) != len(output_frame):
        raise Stage6InputError(
            "Stage 6 changed the number of variant rows: "
            f"{len(input_frame)} input, {len(output_frame)} output."
        )
    expected_columns = list(input_frame.columns) + STAGE6_OUTPUT_COLUMNS
    if list(output_frame.columns) != expected_columns:
        raise Stage6InputError("Stage 6 output columns do not match the approved contract.")

    input_index = build_variant_index(input_frame)
    output_index = build_variant_index(output_frame)
    if set(input_index) != set(output_index):
        raise Stage6InputError("Stage 6 output changed the set of variant keys.")
    input_keys = list(input_index)
    output_order = [output_index[key] for key in input_keys]
    for column in input_frame.columns:
        input_values = input_frame[column].tolist()
        output_values = output_frame[column].iloc[output_order].tolist()
        mismatch = next(
            (
                position
                for position, (input_value, output_value) in enumerate(
                    zip(input_values, output_values, strict=True)
                )
                if str(input_value) != str(output_value)
            ),
            None,
        )
        if mismatch is not None:
            raise Stage6InputError(
                f"Stage 6 changed original field {column!r} "
                f"for {input_keys[mismatch]}."
            )

    for row in output_frame.to_dict(orient="records"):
        expected = calculate_final_score(*[float(row[name]) for name in ["Variant_pathogenicity_score", "Inheritance_score", "Gene_disease_evidence_score", "Case_relevance_score", "Rarity_score"]])
        if abs(float(row["Final_score"]) - expected) > 0.011:
            raise Stage6InputError("Final score disagrees with component weights")
        if to_bool(row.get("AD_fallback")) is True:
            if row.get("Scored_inheritance") != "AD" or not is_missing(row.get("Scored_compound_partner")) or not is_missing(row.get("Compound_pair_score")):
                raise Stage6InputError("AD fallback has inconsistent selected hypothesis fields")
            compound_genes = {gene.casefold() for gene in parse_compound_genes(row.get("Compound_het_gene"))}
            if str(row.get("Scored_gene", "")).casefold() in compound_genes:
                raise Stage6InputError("AD fallback cannot use the same AR compound gene")

    pair_columns = [
        "Chr_Pos",
        "Variant_REF_ALT",
        "Scored_gene",
        "Scored_compound_partner",
        "Compound_pair_score",
    ]
    output_by_key = {
        variant_key(row): row
        for row in output_frame[pair_columns].to_dict(orient="records")
    }
    input_gene_symbols = input_frame["Gene_symbol"].tolist()
    for input_position, key in enumerate(input_keys):
        output_row = output_by_key[key]
        scored_gene = str(output_row["Scored_gene"]).strip()
        row_genes = {
            gene.casefold()
            for gene in parse_gene_symbols(input_gene_symbols[input_position])
        }
        if scored_gene and scored_gene.casefold() not in row_genes:
            raise Stage6InputError(
                f"Scored_gene {scored_gene!r} is not represented by {key}."
            )

    ranges = {
        "Variant_pathogenicity_score": (0.0, 100.0),
        "Pathogenicity_evidence_coverage": (0.0, 1.0),
        "Inheritance_score": (0.0, 100.0),
        "Gene_disease_evidence_score": (0.0, 100.0),
        "Case_relevance_score": (0.0, 100.0),
        "Rarity_score": (0.0, 100.0),
        "Final_score": (0.0, 100.0),
    }
    for column, (minimum, maximum) in ranges.items():
        values = pd.to_numeric(output_frame[column], errors="coerce")
        if values.isna().any() or not values.between(minimum, maximum).all():
            raise Stage6InputError(
                f"Stage 6 output column {column} contains invalid scores."
            )

    present_pair_scores = output_frame.loc[
        ~output_frame["Compound_pair_score"].map(is_missing),
        "Compound_pair_score",
    ]
    numeric_pair_scores = pd.to_numeric(present_pair_scores, errors="coerce")
    if numeric_pair_scores.isna().any() or not numeric_pair_scores.between(
        0.0,
        100.0,
    ).all():
        raise Stage6InputError(
            "Stage 6 output column Compound_pair_score contains invalid scores."
        )

    selected_pairs_by_gene: dict[str, set[frozenset[str]]] = {}
    for key, row in output_by_key.items():
        partner_key = row.get("Scored_compound_partner")
        pair_score = row.get("Compound_pair_score")
        if is_missing(partner_key):
            if not is_missing(pair_score):
                raise Stage6InputError(
                    f"{key} has Compound_pair_score without a selected partner."
                )
            continue
        if is_missing(pair_score):
            raise Stage6InputError(
                f"{key} has a selected compound partner without Compound_pair_score."
            )
        partner_key = str(partner_key).strip()
        if partner_key not in output_by_key:
            raise Stage6InputError(
                f"Selected compound partner {partner_key} for {key} is absent."
            )
        partner = output_by_key[partner_key]
        if str(partner.get("Scored_compound_partner", "")).strip() != key:
            raise Stage6InputError(
                f"Selected compound pair {key} / {partner_key} is not reciprocal."
            )
        if to_float(pair_score) != to_float(partner.get("Compound_pair_score")):
            raise Stage6InputError(
                f"Selected compound pair {key} / {partner_key} has unequal scores."
            )
        if str(row.get("Scored_gene", "")).casefold() != str(
            partner.get("Scored_gene", "")
        ).casefold():
            raise Stage6InputError(
                f"Selected compound pair {key} / {partner_key} has different genes."
            )
        gene_key = str(row.get("Scored_gene", "")).casefold()
        selected_pairs_by_gene.setdefault(gene_key, set()).add(
            frozenset((key, partner_key))
        )

    genes_with_multiple_pairs = [
        gene
        for gene, pairs in selected_pairs_by_gene.items()
        if len(pairs) > 1
    ]
    if genes_with_multiple_pairs:
        raise Stage6InputError(
            "Stage 6 selected more than one compound pair for gene(s): "
            + ", ".join(sorted(genes_with_multiple_pairs))
        )


    mitochondrial_count = sum(is_mitochondrial_variant(row) for row in input_frame.to_dict(orient="records"))
    warnings = ["Partial mitochondrial scoring: MT-specific pathogenicity and rarity unavailable; existing missing rarity default 50/100."] if mitochondrial_count else []
    record_qc("stage6_score", len(input_frame), len(output_frame), len(input_index), len(output_index), "rank", warnings=warnings, partially_scored_mitochondrial_variants=mitochondrial_count)

def print_scoring_summary(
    input_frame: pd.DataFrame,
    output_frame: pd.DataFrame,
    warnings: Sequence[str],
    selected_candidates: Sequence[Candidate],
) -> None:
    genes_per_row = input_frame["Gene_symbol"].map(
        lambda value: len(parse_gene_symbols(value))
    )
    transcripts_per_row = input_frame["Transcript_ID"].map(
        lambda value: len(split_aligned(value))
    )
    print(f"Input variants: {len(input_frame):,}")
    print(f"Output variants: {len(output_frame):,}")
    print(f"Unique scored genes: {output_frame['Scored_gene'].replace('', pd.NA).nunique():,}")
    print(f"Rows with one gene: {(genes_per_row == 1).sum():,}")
    print(f"Rows with multiple genes: {(genes_per_row > 1).sum():,}")
    print(f"Rows with one transcript: {(transcripts_per_row == 1).sum():,}")
    print(f"Rows with multiple transcripts: {(transcripts_per_row > 1).sum():,}")
    print(f"Case-gene winners: {sum(item.is_case_gene for item in selected_candidates):,}")
    print(f"Missing Scored_gene: {output_frame['Scored_gene'].map(is_missing).sum():,}")
    print(f"Warnings / parsing failures: {len(warnings):,}")
    print(
        "Retained compound-het alleles: "
        f"{(~output_frame['Scored_compound_partner'].map(is_missing)).sum():,}"
    )

    score_columns = [
        "Variant_pathogenicity_score",
        "Pathogenicity_evidence_coverage",
        "Inheritance_score",
        "Gene_disease_evidence_score",
        "Case_relevance_score",
        "Rarity_score",
        "Compound_pair_score",
        "Final_score",
    ]
    print("\nScore distributions:")
    print(describe_scores(output_frame, score_columns).round(2).to_string())

    final_values = pd.to_numeric(output_frame["Final_score"], errors="coerce")
    histogram = pd.cut(
        final_values,
        bins=[-0.001, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100.001],
        right=False,
        labels=[
            "0-9",
            "10-19",
            "20-29",
            "30-39",
            "40-49",
            "50-59",
            "60-69",
            "70-79",
            "80-89",
            "90-100",
        ],
    ).value_counts(sort=False)
    print("\nFinal_score histogram:")
    print(histogram.to_string())

    top_columns = [
        "Chr_Pos",
        "Variant_REF_ALT",
        "Scored_gene",
        "Scored_compound_partner",
        "Compound_pair_score",
        "Final_score",
        "Variant_pathogenicity_score",
        "Inheritance_score",
        "Gene_disease_evidence_score",
        "Case_relevance_score",
        "Rarity_score",
        "Case_gene_match",
    ]
    print("\nTop candidates:")
    print(output_frame[top_columns].head(20).to_string(index=False))


def print_validation_summary(
    frame: pd.DataFrame,
    warnings: Sequence[str],
    gene_cache: Mapping[str, Any],
    variant_index: Mapping[str, int],
) -> None:
    print(f"Input variants: {len(frame):,}")
    print(f"Input columns: {len(frame.columns):,}")
    print(f"Cached genes: {len(gene_cache):,}")
    print(f"Indexed variants: {len(variant_index):,}")
    print(f"Parsing warnings: {len(warnings):,}")
    for warning in warnings[:10]:
        print(f"  Warning: {warning}")
    if len(warnings) > 10:
        print(f"  ... {len(warnings) - 10:,} additional warnings")
    print("Structural validation passed.")


def run(args: argparse.Namespace) -> None:
    frame = load_stage5(args.input)
    context = load_case_context(args.case_context)
    warnings = validate_inputs(frame, context)
    gene_cache = build_gene_cache(frame, context)
    variant_index = build_variant_index(frame)
    print_validation_summary(frame, warnings, gene_cache, variant_index)

    if args.validate_only:
        return
    if args.output is None:
        raise Stage6InputError("--output is required unless --validate-only is used.")

    output, selected = score_stage6(
        frame,
        context,
        gene_cache,
        variant_index,
        warnings,
        progress_every=args.progress_every,
    )
    validate_scored_output(frame, output)
    print("Writing...")
    atomic_write_csv(output, args.output)
    print_scoring_summary(frame, output, warnings, selected)


if __name__ == "__main__":
    if "snakemake" in globals():
        run(args_from_snakemake(snakemake))
    else:
        run(parse_args())
