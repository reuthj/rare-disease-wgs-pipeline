#!/usr/bin/env python3
"""Add local Stage 5 gene evidence without changing variant rows."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import tempfile
from collections import defaultdict
from pathlib import Path

import pandas as pd
from workflow_checks import validate_table, check_transition
import yaml


LOCAL_COLUMNS = [
    "Disease_associations",
    "Missense_constraint",
    "LOEUF",
    "ClinVar_pathogenic_fraction",
    "ClinVar_pathogenic_count",
    "Known_variant_landscape",
    "Gene_HPO_terms",
    "Matched_HPO_terms",
    "Inheritance_and_mechanism",
]

GENE_COLUMN_ALIASES = ("gene_symbol", "Gene", "SYMBOL", "gene")
ALIGNED_GENE_COLUMN = "Gene_symbol_transcript_aligned"
CLINVAR_CACHE_VERSION = 2


def _find_column(columns, aliases, source, required=True):
    """Return the first matching column, allowing harmless case differences."""
    exact = {str(column): column for column in columns}
    lower = {str(column).lower(): column for column in columns}
    for alias in aliases:
        if alias in exact:
            return exact[alias]
        if alias.lower() in lower:
            return lower[alias.lower()]
    if required:
        raise ValueError(f"{source} is missing a required column: {aliases}")
    return None


def parse_gene_symbols(value):
    """Return distinct genes in transcript order, skipping empty positions."""
    parts = [
        re.sub(r"^-\s+(?=\S)", "", part.strip()).strip().upper()
        for part in str(value).split("|")
    ]
    return list(dict.fromkeys(
        part for part in parts if part and part not in {"-", "."}
    ))


def clean_gene_symbol(value):
    return " | ".join(parse_gene_symbols(value))


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Stage 4 CSV")
    parser.add_argument("--output", required=True, help="Stage 5 local CSV")
    parser.add_argument("--case-context", required=True)
    parser.add_argument("--gencc", required=True)
    parser.add_argument("--constraint", required=True)
    parser.add_argument("--clinvar", required=True)
    parser.add_argument("--hpo-annotations", required=True)
    parser.add_argument("--hpo-ontology", required=True)
    parser.add_argument(
        "--cache-dir",
        default="results/intermediate/stage5/cache",
        help="Small reusable source summaries; use an empty value to disable",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=500_000,
        help="Report long scans every N records; 0 disables counters",
    )
    return parser.parse_args(argv)


def load_gencc_evidence(path, genes):
    """Return source-specific disease assertions and explicit GenCC MOI."""
    result = {gene: {} for gene in genes}
    if not genes:
        return result

    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    gene_col = _find_column(
        table.columns,
        ("gene_symbol", "gene", "submitted_as_hgnc_symbol"),
        "GenCC",
    )
    disease_col = _find_column(
        table.columns,
        ("disease_title", "disease_original_title", "submitted_as_disease_name"),
        "GenCC",
    )
    disease_id_col = _find_column(
        table.columns,
        ("disease_curie", "disease_original_curie", "submitted_as_disease_id"),
        "GenCC",
        required=False,
    )
    class_col = _find_column(
        table.columns,
        ("classification_title", "submitted_as_classification_name"),
        "GenCC",
    )
    submitter_col = _find_column(
        table.columns,
        ("submitter_title", "submitted_as_submitter_name"),
        "GenCC",
    )
    moi_col = _find_column(
        table.columns,
        ("moi_title", "submitted_as_moi_name"),
        "GenCC",
        required=False,
    )

    table[gene_col] = table[gene_col].str.strip().str.upper()
    table = table[table[gene_col].isin(genes)]
    # Orphanet is not a Stage 5 evidence source, even when bundled in GenCC.
    table = table[
        table[submitter_col].str.strip().str.casefold() != "orphanet"
    ]

    for gene, group in table.groupby(gene_col, sort=False):
        assertions = defaultdict(set)
        inheritance = defaultdict(set)
        diseases = set()

        for row in group.to_dict("records"):
            disease = row[disease_col].strip()
            classification = row[class_col].strip()
            submitter = row[submitter_col].strip()
            disease_id = row.get(disease_id_col, "").strip() if disease_id_col else ""
            moi = row.get(moi_col, "").strip() if moi_col else ""
            if not disease or not classification or not submitter:
                continue

            assertions[disease].add(classification)
            diseases.add((disease_id, disease))
            if moi:
                inheritance[disease].add(moi)

        disease_parts = []
        classification_order = {
            "definitive": 0,
            "strong": 1,
            "moderate": 2,
            "supportive": 3,
            "limited": 4,
            "disputed evidence": 5,
            "disputed": 5,
            "refuted evidence": 6,
            "refuted": 6,
            "no known disease relationship": 7,
            "animal model only": 8,
        }
        for disease in sorted(assertions, key=str.casefold):
            classification = min(
                assertions[disease],
                key=lambda value: (
                    classification_order.get(value.casefold(), 99),
                    value.casefold(),
                ),
            )
            disease_parts.append(f"{disease} - {classification}")

        inheritance_parts = [
            f"{disease}: {'; '.join(sorted(inheritance[disease]))}"
            for disease in sorted(inheritance, key=str.casefold)
        ]
        result[gene] = {
            "Disease_associations": " || ".join(disease_parts),
            "Inheritance_and_mechanism": " || ".join(inheritance_parts),
            "_diseases": sorted(diseases, key=lambda item: (item[1].casefold(), item[0])),
        }

    return result


def load_constraint_evidence(path, genes, progress_every=0):
    """Stream gnomAD and retain only the best usable row per candidate gene."""
    result = {
        gene: {
            "Missense_constraint": "",
            "LOEUF": "",
        }
        for gene in genes
    }
    if not genes:
        return result

    compression = "gzip" if str(path).endswith((".gz", ".bgz")) else "infer"
    header = pd.read_csv(path, sep="\t", nrows=0, compression=compression)
    gene_col = _find_column(header.columns, ("gene", "gene_symbol"), "gnomAD")
    transcript_col = _find_column(header.columns, ("transcript",), "gnomAD")
    oe_col = _find_column(header.columns, ("mis.oe", "mis_oe"), "gnomAD")
    z_col = _find_column(
        header.columns,
        ("mis.z_score", "mis_z_score", "mis.z"),
        "gnomAD",
    )
    loeuf_col = _find_column(
        header.columns,
        (
            "lof.oe_ci.upper",
            "lof_oe_ci_upper",
            "lof.oe_ci_upper",
        ),
        "gnomAD",
    )
    mane_col = _find_column(
        header.columns, ("mane_select", "mane"), "gnomAD", required=False
    )
    canonical_col = _find_column(
        header.columns, ("canonical",), "gnomAD", required=False
    )
    flags_col = _find_column(
        header.columns, ("constraint_flags",), "gnomAD", required=False
    )

    def preferred(value):
        value = str(value).strip().casefold()
        return value not in {"", "0", "false", "no", "n", "na", "nan", "none"}

    usecols = [
        column
        for column in (
            gene_col,
            transcript_col,
            oe_col,
            z_col,
            loeuf_col,
            mane_col,
            canonical_col,
            flags_col,
        )
        if column
    ]
    best = {}
    rows_read = 0
    next_report = progress_every
    for chunk in pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        compression=compression,
        usecols=list(dict.fromkeys(usecols)),
        chunksize=100_000,
    ):
        rows_read += len(chunk)
        chunk[gene_col] = chunk[gene_col].str.strip().str.upper()
        chunk = chunk[chunk[gene_col].isin(genes)]
        for row in chunk.to_dict("records"):
            try:
                oe = float(row[oe_col])
                z_score = float(row[z_col])
                loeuf = float(row[loeuf_col])
            except (TypeError, ValueError):
                continue
            rank = (
                0 if mane_col and preferred(row[mane_col]) else 1,
                0 if canonical_col and preferred(row[canonical_col]) else 1,
                str(row[transcript_col]),
            )
            gene = row[gene_col]
            if gene not in best or rank < best[gene][0]:
                best[gene] = (rank, oe, z_score, loeuf, row.get(flags_col, ""))

        if progress_every and rows_read >= next_report:
            print(f"  gnomAD rows scanned: {rows_read:,}")
            next_report += progress_every

    for gene, (_, oe, z_score, loeuf, flags) in best.items():
        missense_value = f"O/E={oe:.3g} | Z={z_score:.3g}"

        flags = str(flags).strip()
        if flags and flags.casefold() not in {"[]", "na", "nan", "none"}:
            missense_value += f" | flag={flags}"

        result[gene]["Missense_constraint"] = missense_value
        result[gene]["LOEUF"] = f"{loeuf:.3g}"

    if progress_every:
        print(f"  gnomAD scan complete: {rows_read:,} rows")

    return result


def load_clinvar_evidence(path, genes, cache_dir=None, progress_every=0):
    """Stream ClinVar once; reuse a small summary when source and genes match."""
    empty = {
        "ClinVar_pathogenic_fraction": "",
        "ClinVar_pathogenic_count": 0,
        "Known_variant_landscape": "",
    }
    result = {gene: empty.copy() for gene in genes}
    if not genes:
        return result

    source_path = Path(path)
    stat = source_path.stat()
    cache_path = None
    if cache_dir:
        cache_key = hashlib.sha256(
            json.dumps(
                {
                    "version": CLINVAR_CACHE_VERSION,
                    "source": source_path.name,
                    "size": stat.st_size,
                    "mtime_ns": stat.st_mtime_ns,
                    "genes": sorted(genes),
                },
                sort_keys=True,
            ).encode()
        ).hexdigest()[:20]
        cache_path = Path(cache_dir) / f"clinvar_{cache_key}.json"
        if cache_path.is_file():
            try:
                with open(cache_path, encoding="utf-8") as handle:
                    cached = json.load(handle)
                if cached.get("version") == CLINVAR_CACHE_VERSION:
                    cached_evidence = cached.get("evidence", {})
                    if set(cached_evidence) == genes:
                        print(f"  ClinVar cache hit: {cache_path.name}")
                        return cached_evidence
            except (OSError, json.JSONDecodeError, TypeError) as error:
                print(f"WARNING - ignoring unusable ClinVar cache: {error}")

    pathogenic_terms = {
        "pathogenic",
        "likely pathogenic",
        "pathogenic/likely pathogenic",
        "likely pathogenic/pathogenic",
        "pathogenic, low penetrance",
        "likely pathogenic, low penetrance",
    }
    eligible_terms = pathogenic_terms | {
        "benign",
        "likely benign",
        "benign/likely benign",
        "likely benign/benign",
        "uncertain significance",
        "vus-high",
        "vus-mid",
        "vus-low",
        "vus high",
        "vus mid",
        "vus low",
    }
    category_rules = (
        ("splice", ("splice_acceptor_variant", "splice_donor_variant", "splice_region_variant")),
        ("frameshift", ("frameshift_variant",)),
        ("nonsense", ("stop_gained", "nonsense")),
        ("stop-loss", ("stop_lost",)),
        ("start-loss", ("start_lost",)),
        ("missense", ("missense_variant",)),
        ("in-frame", ("inframe_insertion", "inframe_deletion")),
        ("synonymous", ("synonymous_variant",)),
        ("UTR", ("utr_variant",)),
        ("intronic", ("intron_variant",)),
        ("noncoding RNA", ("non_coding_transcript", "ncrna")),
        ("near-gene", ("upstream_gene_variant", "downstream_gene_variant", "intergenic_variant")),
        ("regulatory", ("regulatory_region", "tf_binding_site")),
        ("transcript change", ("transcript_ablation", "transcript_amplification")),
    )
    category_rank = {name: index for index, (name, _) in enumerate(category_rules)}
    variants = {gene: {} for gene in genes}
    records_scanned = 0
    candidate_records = 0
    next_report = progress_every

    opener = gzip.open if str(path).endswith((".gz", ".bgz")) else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            records_scanned += 1
            if progress_every and records_scanned >= next_report:
                print(
                    f"  ClinVar records scanned: {records_scanned:,}; "
                    f"candidate records: {candidate_records:,}"
                )
                next_report += progress_every
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8 or fields[2] in {"", "."}:
                continue

            info = {}
            for item in fields[7].split(";"):
                key, separator, value = item.partition("=")
                if separator:
                    info[key] = value

            gene_info = info.get("GENEINFO", "")
            record_genes = {
                token.split(":", 1)[0].strip().upper()
                for token in re.split(r"[|,]", gene_info)
                if token.strip()
            } & genes
            if not record_genes:
                continue
            candidate_records += 1

            clinsig = re.sub(r"\s+", " ", info.get("CLNSIG", "").replace("_", " "))
            clinsig = clinsig.strip().casefold()
            if clinsig not in eligible_terms:
                continue
            status = "pathogenic" if clinsig in pathogenic_terms else "eligible"

            consequence_text = info.get("MC", "").casefold()
            category = ""
            for name, terms in category_rules:
                if any(term in consequence_text for term in terms):
                    category = name
                    break
            if not category:
                raw_terms = [
                    item.split("|", 1)[-1].strip().replace("_", " ")
                    for item in consequence_text.split(",")
                    if item.strip()
                ]
                category = raw_terms[0] if raw_terms else "unclassified"

            variation_id = fields[2]
            for gene in record_genes:
                previous = variants[gene].get(variation_id)
                if previous is None:
                    variants[gene][variation_id] = [status, category]
                    continue
                if status == "pathogenic":
                    previous[0] = "pathogenic"
                if category_rank.get(category, 99) < category_rank.get(previous[1], 99):
                    previous[1] = category

    for gene, gene_variants in variants.items():
        eligible_count = len(gene_variants)
        pathogenic = [
            values for values in gene_variants.values() if values[0] == "pathogenic"
        ]
        pathogenic_count = len(pathogenic)
        result[gene]["ClinVar_pathogenic_count"] = pathogenic_count
        if eligible_count:
            result[gene]["ClinVar_pathogenic_fraction"] = round(
                pathogenic_count / eligible_count, 6
            )

        landscape_counts = defaultdict(int)
        for _, category in pathogenic:
            landscape_counts[category] += 1
        if landscape_counts:
            known = [name for name, _ in category_rules]
            ordered = known + sorted(set(landscape_counts) - set(known))
            parts = [
                f"{name}={landscape_counts[name]}"
                for name in ordered
                if landscape_counts[name]
            ]
            result[gene]["Known_variant_landscape"] = "P/LP: " + " | ".join(parts)

    if progress_every:
        print(
            f"  ClinVar scan complete: {records_scanned:,} records; "
            f"candidate records: {candidate_records:,}"
        )
    if cache_path:
        temporary_path = None
        try:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(
                mode="w",
                suffix=".tmp",
                dir=cache_path.parent,
                delete=False,
                encoding="utf-8",
            ) as temporary:
                json.dump(
                    {"version": CLINVAR_CACHE_VERSION, "evidence": result},
                    temporary,
                    sort_keys=True,
                )
                temporary_path = Path(temporary.name)
            os.replace(temporary_path, cache_path)
        except OSError as error:
            print(f"WARNING - ClinVar cache was not saved: {error}")
        finally:
            if temporary_path and temporary_path.exists():
                temporary_path.unlink()

    return result


def load_hpo_evidence(
    annotation_path,
    ontology_path,
    genes,
    gencc,
    case_context,
    progress_every=0,
):
    """Build disease HPO profiles and exact/ancestor patient matches."""
    empty = {"Gene_HPO_terms": "", "Matched_HPO_terms": ""}
    result = {gene: empty.copy() for gene in genes}
    if not genes:
        return result

    def normalize_hp(value):
        match = re.search(r"HP[_:](\d{7})", str(value).upper())
        return f"HP:{match.group(1)}" if match else ""

    def normalize_disease_id(value):
        value = str(value).strip().upper()
        if value.startswith("MIM:"):
            value = "OMIM:" + value.removeprefix("MIM:")
        return value

    def normalize_name(value):
        return re.sub(r"[^a-z0-9]+", " ", str(value).casefold()).strip()

    def context_values(key):
        values = case_context.get(key) or []
        return [values] if isinstance(values, str) else values

    with open(ontology_path, encoding="utf-8") as handle:
        ontology = json.load(handle)
    graphs = ontology.get("graphs") or []
    if not graphs:
        raise ValueError("HPO ontology JSON does not contain graphs")

    labels = {}
    parents = defaultdict(set)
    for graph in graphs:
        for node in graph.get("nodes", []):
            term = normalize_hp(node.get("id", ""))
            if term:
                labels[term] = str(node.get("lbl", "")).strip()
        for edge in graph.get("edges", []):
            child = normalize_hp(edge.get("sub", ""))
            parent = normalize_hp(edge.get("obj", ""))
            predicate = str(edge.get("pred", "")).casefold()
            if child and parent and (predicate == "is_a" or "subclassof" in predicate):
                parents[child].add(parent)

    genes_by_id = defaultdict(set)
    genes_by_name = defaultdict(set)
    case_genes = {
        str(gene).strip().upper() for gene in context_values("case_genes")
    }
    case_ids = {
        normalize_disease_id(value)
        for key in ("orpha_ids", "mondo_ids", "disease_ids")
        for value in context_values(key)
        if normalize_disease_id(value)
    }

    for gene in genes:
        diseases = set(gencc.get(gene, {}).get("_diseases", []))
        if gene in case_genes:
            diseases.update((disease_id, "") for disease_id in case_ids)
        for disease_id, disease_name in diseases:
            disease_id = normalize_disease_id(disease_id)
            disease_name = normalize_name(disease_name)
            if disease_id:
                genes_by_id[disease_id].add(gene)
            if disease_name:
                genes_by_name[disease_name].add(gene)

    profiles = {gene: defaultdict(set) for gene in genes}
    with open(annotation_path, encoding="utf-8") as handle:
        lines = (
            line[1:] if line.lower().startswith("#database_id") else line
            for line in handle
            if line.strip()
            and (not line.startswith("#") or line.lower().startswith("#database_id"))
        )
        reader = csv.DictReader(lines, delimiter="\t")
        fields = reader.fieldnames or []
        disease_id_col = _find_column(
            fields, ("database_id", "disease_id"), "HPO annotations"
        )
        disease_name_col = _find_column(
            fields, ("disease_name",), "HPO annotations"
        )
        qualifier_col = _find_column(
            fields, ("qualifier",), "HPO annotations", required=False
        )
        hpo_col = _find_column(fields, ("hpo_id",), "HPO annotations")
        aspect_col = _find_column(
            fields, ("aspect",), "HPO annotations", required=False
        )

        rows_read = 0
        next_report = progress_every
        for row in reader:
            rows_read += 1
            if progress_every and rows_read >= next_report:
                print(f"  HPO annotation rows scanned: {rows_read:,}")
                next_report += progress_every
            if qualifier_col and row.get(qualifier_col, "").strip().upper() == "NOT":
                continue
            if aspect_col and row.get(aspect_col, "").strip().upper() != "P":
                continue

            disease_id = normalize_disease_id(row.get(disease_id_col, ""))
            disease_name = row.get(disease_name_col, "").strip()
            hpo_id = normalize_hp(row.get(hpo_col, ""))
            if not hpo_id:
                continue

            matched_genes = genes_by_id.get(disease_id, set())
            if not matched_genes:
                matched_genes = genes_by_name.get(normalize_name(disease_name), set())
            for gene in matched_genes:
                profiles[gene][(disease_id, disease_name)].add(hpo_id)

        if progress_every:
            print(f"  HPO scan complete: {rows_read:,} annotation rows")

    patient_terms = {
        normalize_hp(value)
        for value in context_values("hpo_terms")
        if normalize_hp(value)
    }
    ancestor_cache = {}

    def ancestors(term):
        if term not in ancestor_cache:
            found = set()
            stack = list(parents.get(term, set()))
            while stack:
                parent = stack.pop()
                if parent in found:
                    continue
                found.add(parent)
                stack.extend(parents.get(parent, set()))
            ancestor_cache[term] = found
        return ancestor_cache[term]

    def term_text(term):
        return f"{term} {labels.get(term, '')}".strip()

    for gene, diseases in profiles.items():
        groups = []
        disease_terms = set()
        for (disease_id, disease_name), terms in sorted(
            diseases.items(), key=lambda item: (item[0][1].casefold(), item[0][0])
        ):
            disease_terms.update(terms)
            title = disease_name or disease_id
            if disease_id and disease_name:
                title += f" ({disease_id})"
            groups.append(
                f"{title}: " + "; ".join(term_text(term) for term in sorted(terms))
            )
        result[gene]["Gene_HPO_terms"] = " || ".join(groups)

        matches = set()
        for patient_term in patient_terms:
            for disease_term in disease_terms:
                if (
                    patient_term == disease_term
                    or disease_term in ancestors(patient_term)
                    or patient_term in ancestors(disease_term)
                ):
                    matches.add(
                        f"{term_text(patient_term)} ↔ {term_text(disease_term)}"
                    )
        result[gene]["Matched_HPO_terms"] = " | ".join(sorted(matches))

    return result


def main(argv=None):
    args = parse_args(argv)
    input_path = Path(args.input)
    output_path = Path(args.output)
    context_path = Path(args.case_context)
    if not input_path.is_file():
        raise FileNotFoundError(f"Stage 4 input not found: {input_path}")
    if not context_path.is_file():
        raise FileNotFoundError(f"Case context not found: {context_path}")

    print("Stage 5 local evidence: loading input")
    input_df = pd.read_csv(input_path, dtype=str, keep_default_na=False)
    validate_table(input_df, "stage5_gene_evidence input", required=["MANE_Transcript_ID"], aligned=True)
    original_df = input_df.copy(deep=True)
    overlap = [
        column for column in [ALIGNED_GENE_COLUMN, *LOCAL_COLUMNS]
        if column in input_df.columns
    ]
    if overlap:
        raise ValueError(
            "Stage 5 columns already exist in the input: " + ", ".join(overlap)
        )

    gene_col = _find_column(
        input_df.columns, GENE_COLUMN_ALIASES, "Stage 4 input"
    )
    gene_lists = input_df[gene_col].map(parse_gene_symbols)
    gene_keys = gene_lists.map(" | ".join)
    changed = (gene_keys != input_df[gene_col]).sum()
    print(f"  Gene symbols cleaned: {changed:,}/{len(input_df):,}")
    print(f"  Rows with multiple genes: {(gene_lists.map(len) > 1).sum():,}")
    genes = {gene for row_genes in gene_lists for gene in row_genes}
    with open(context_path, encoding="utf-8") as handle:
        case_context = yaml.safe_load(handle) or {}
    if not isinstance(case_context, dict):
        raise ValueError("Case context must be a YAML mapping")
    progress_every = max(0, args.progress_every)
    cache_dir = args.cache_dir or None

    evidence = {
        gene: {column: "" for column in LOCAL_COLUMNS}
        for gene in genes
    }
    source_errors = []

    def run_source(label, function, *function_args):
        print(f"Stage 5 local evidence: {label}")
        try:
            return function(*function_args)
        except Exception as error:  # keep independent enrichment failures isolated
            source_errors.append(f"{label}: {error}")
            print(f"ERROR - {label} unavailable; its fields will remain blank: {error}")
            return {}

    gencc = run_source("GenCC lookup", load_gencc_evidence, args.gencc, genes)
    constraint = run_source(
        "gnomAD constraint lookup",
        load_constraint_evidence,
        args.constraint,
        genes,
        progress_every,
    )
    clinvar = run_source(
        "ClinVar scan",
        load_clinvar_evidence,
        args.clinvar,
        genes,
        cache_dir,
        progress_every,
    )
    hpo = run_source(
        "HPO matching",
        load_hpo_evidence,
        args.hpo_annotations,
        args.hpo_ontology,
        genes,
        gencc,
        case_context,
        progress_every,
    )

    for gene in genes:
        for source in (gencc, constraint, clinvar, hpo):
            values = source.get(gene, {})
            for column in LOCAL_COLUMNS:
                if column in values:
                    evidence[gene][column] = values[column]

    print("Stage 5 local evidence: assembling output")
    output_df = input_df.copy(deep=True)
    output_df[gene_col] = gene_keys
    output_df.insert(
        input_df.columns.get_loc(gene_col) + 1,
        ALIGNED_GENE_COLUMN,
        original_df[gene_col],
    )
    for column in LOCAL_COLUMNS:
        values = []
        for row_genes in gene_lists:
            if len(row_genes) > 1:
                value = " || ".join(
                    f"{gene}: {evidence[gene][column]}"
                    for gene in row_genes
                )
            elif row_genes:
                value = evidence[row_genes[0]][column]
            else:
                value = ""
            values.append(value)
        output_df[column] = values

    if len(output_df) != len(input_df):
        raise AssertionError("Stage 5 changed the number of variant rows")
    if not output_df[input_df.columns.drop(gene_col)].equals(
        original_df.drop(columns=[gene_col])
    ):
        raise AssertionError("Stage 5 changed a pre-existing value outside the gene column")
    if not output_df[gene_col].equals(gene_keys):
        raise AssertionError("Stage 5 gene column does not match cleaned symbols")
    if not output_df[ALIGNED_GENE_COLUMN].equals(original_df[gene_col]):
        raise AssertionError("Stage 5 changed the original transcript-aligned genes")
    expected_columns = list(input_df.columns)
    expected_columns.insert(input_df.columns.get_loc(gene_col) + 1, ALIGNED_GENE_COLUMN)
    if list(output_df.columns) != expected_columns + LOCAL_COLUMNS:
        raise AssertionError("Unexpected Stage 5 output columns")

    check_transition(original_df, output_df, "stage5_gene_evidence", aligned=True, changed_fields=(gene_col,), warnings=source_errors)

    print("Stage 5 local evidence: writing output")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".tmp",
        prefix=f".{output_path.name}.",
        dir=output_path.parent,
        delete=False,
    ) as temporary:
        temporary_path = Path(temporary.name)
    try:
        output_df.to_csv(temporary_path, index=False, encoding="utf-8-sig")
        os.replace(temporary_path, output_path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()

    matched = {
        "GenCC": sum(bool(evidence[g]["Disease_associations"]) for g in genes),
        "gnomAD constraint": sum(bool(evidence[g]["Missense_constraint"]) for g in genes),
        "ClinVar classified": sum(
            evidence[g]["ClinVar_pathogenic_fraction"] != "" for g in genes
        ),
        "HPO profiles": sum(bool(evidence[g]["Gene_HPO_terms"]) for g in genes),
    }
    print("\nStage 5 local evidence summary")
    print(f"Input variants:  {len(input_df):,}")
    print(f"Output variants: {len(output_df):,}")
    print(f"Unique genes:    {len(genes):,}")
    for label, count in matched.items():
        print(f"{label}: {count:,}/{len(genes):,}")
    print(f"Source errors: {len(source_errors)}")


if __name__ == "__main__":
    main()
