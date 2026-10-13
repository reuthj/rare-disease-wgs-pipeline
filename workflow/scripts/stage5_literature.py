#!/usr/bin/env python3
"""Add Europe PMC literature evidence without changing existing variant data."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import tempfile
import time
from pathlib import Path

import pandas as pd
from workflow_checks import validate_table, check_transition
import requests
import yaml


API_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
CACHE_VERSION = 4
PAGE_SIZE = 5
MAX_REVIEWED = 8
MAX_TOP_PUBLICATIONS = 3

OUTPUT_COLUMNS = ["Literature_support", "Top_publications"]
GENE_COLUMN_ALIASES = ("Gene_symbol", "gene_symbol", "Gene", "SYMBOL", "gene")
ALIGNED_GENE_COLUMN = "Gene_symbol_transcript_aligned"
HGVSC_ALIASES = ("HGVSc", "HGVS_c", "hgvsc")
HGVSP_ALIASES = ("HGVSp", "HGVS_p", "hgvsp")
RSID_ALIASES = ("Existing_variation", "rsID", "Known_variant", "RSID")

STUDY_PATTERNS = {
    "human-cases": re.compile(
        r"\b(patient|patients|proband|probands|family|families|cohort|"
        r"case report|case series|clinical features?|affected individuals?)\b",
        re.IGNORECASE,
    ),
    "functional": re.compile(
        r"\b(functional (?:analysis|assay|study)|assay|enzyme activity|"
        r"electrophysiolog\w*|in vitro|cell(?:ular)? model|protein expression|"
        r"transcript expression|western blot)\b",
        re.IGNORECASE,
    ),
    "animal": re.compile(
        r"\b(mouse|mice|murine|zebrafish|rat|rats|drosophila|animal model)\b",
        re.IGNORECASE,
    ),
}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Stage 5A CSV")
    parser.add_argument("--output", required=True, help="Stage 5 final CSV")
    parser.add_argument("--case-context", required=True)
    parser.add_argument(
        "--cache-dir",
        default="results/intermediate/stage5/cache/europe_pmc",
    )
    parser.add_argument("--progress-every", type=int, default=25)
    parser.add_argument("--max-retries", type=int, default=3, help="Retries after the initial request")
    args = parser.parse_args(argv)
    if args.max_retries < 0:
        parser.error("--max-retries must be nonnegative")
    return args


def load_literature_context(path):
    """Use specific disease names, falling back to broad disease domains."""
    with open(path, encoding="utf-8") as handle:
        context = yaml.safe_load(handle) or {}
    if not isinstance(context, dict):
        raise ValueError("Case context must be a YAML mapping")

    def clean_list(value):
        if isinstance(value, str):
            value = [value]
        if not isinstance(value, list):
            return []
        cleaned = [str(item).strip() for item in value]
        return list(dict.fromkeys(item for item in cleaned if item))

    disease_names = clean_list(context.get("disease_names", []))
    broad_domains = clean_list(context.get("broad_disease_domains", []))
    manual_genes = {
        gene.upper() for gene in clean_list(context.get("manual_genes", []))
    }
    return {
        "disease_terms": disease_names or broad_domains,
        "context_source": "disease_names" if disease_names else "broad_disease_domains",
        "manual_genes": manual_genes,
    }


def extract_gene_variants(row, gene_col, aligned_col, hgvsc_col, hgvsp_col, rsid_col):
    """Return ordered genes and their transcript-aligned exact identifiers."""
    missing = {"", "-", ".", "nan", "none"}

    def split_aligned(value):
        return [part.strip() for part in str(value).split("|")]

    aligned_genes = split_aligned(row.get(aligned_col, row.get(gene_col, "")))
    aligned_genes = [
        re.sub(r"^-\s+(?=\S)", "", gene).strip().upper()
        for gene in aligned_genes
    ]
    row_genes = [
        gene.strip().upper()
        for gene in str(row.get(gene_col, "")).split("|")
        if gene.strip().lower() not in missing
    ]
    genes = list(dict.fromkeys(row_genes or [
        gene for gene in aligned_genes if gene.lower() not in missing
    ]))
    identifiers = {gene: [] for gene in genes}

    def add(gene, value, kind):
        if gene not in identifiers or str(value).strip().lower() in missing:
            return
        if kind == "rsid":
            terms = re.findall(r"\brs\d+\b", str(value), flags=re.IGNORECASE)
        else:
            terms = []
            for raw in re.split(r"[&,;]", str(value)):
                raw = raw.strip()
                if not raw or raw.lower() in missing:
                    continue
                terms.append(raw)
                suffix = re.search(r"(?:^|:)((?:c|p)\.[^\s]+)$", raw, re.IGNORECASE)
                if suffix:
                    terms.append(suffix.group(1))
        identifiers[gene].extend(terms)

    for column, kind in ((hgvsc_col, "hgvs"), (hgvsp_col, "hgvs")):
        if not column:
            continue
        parts = split_aligned(row.get(column, ""))
        if len(parts) == len(aligned_genes):
            for gene, value in zip(aligned_genes, parts):
                add(gene, value, kind)
        elif len(genes) == 1:
            add(genes[0], row.get(column, ""), kind)

    if rsid_col:
        rsids = re.findall(
            r"\brs\d+\b", str(row.get(rsid_col, "")), flags=re.IGNORECASE
        )
        for gene in genes:
            for rsid in rsids:
                add(gene, rsid, "rsid")

    for gene in identifiers:
        identifiers[gene] = list(dict.fromkeys(identifiers[gene]))
    return genes, identifiers


def extract_gene_diseases(value, genes):
    """Return Moderate-or-stronger disease associations for each row gene."""
    result = {gene: [] for gene in genes}
    current_gene = genes[0] if len(genes) == 1 else None

    for part in str(value).split(" || "):
        part = part.strip()
        for gene in genes:
            prefix = f"{gene}:"
            if part.upper().startswith(prefix.upper()):
                current_gene = gene
                part = part[len(prefix):].strip()
                break
        if not current_gene or not part:
            continue
        match = re.match(
            r"^(.+?)\s+(?:—|-)\s+(Moderate|Strong|Definitive)\s*$",
            part,
            flags=re.IGNORECASE,
        )
        if match:
            result[current_gene].append(match.group(1).strip())

    for gene in result:
        result[gene] = list(dict.fromkeys(result[gene]))
    return result


def disease_aliases(term):
    """Return a small deterministic set of exact disease-name alternatives."""
    term = re.sub(r"\s+", " ", str(term)).strip(" ,;:-")
    if not term:
        return []

    aliases = []

    def add(value):
        value = re.sub(r"\s+", " ", value).strip(" ,;:-")
        if value and value.casefold() not in {item.casefold() for item in aliases}:
            aliases.append(value)

    add(term)
    without_inheritance = re.sub(
        r"^(?:(?:autosomal|x[- ]linked)\s+(?:dominant|recessive)|"
        r"mitochondrial)\s+",
        "",
        term,
        flags=re.IGNORECASE,
    )
    add(without_inheritance)

    # Add a version without a standalone terminal number.
    # Example: "Hypotrichosis 1" -> "Hypotrichosis".
    # Do not shorten names such as "Bartter disease type 3".
    for value in list(aliases):
        terminal_number = re.search(
            r"\s+\d+[A-Z]?$",
            value,
            flags=re.IGNORECASE,
        )
        if terminal_number:
            shortened = value[:terminal_number.start()].rstrip(" ,;:-")
            if not re.search(r"\btype$", shortened, flags=re.IGNORECASE):
                add(shortened)

    roman = {
        1: "I", 2: "II", 3: "III", 4: "IV", 5: "V",
        6: "VI", 7: "VII", 8: "VIII", 9: "IX", 10: "X",
        11: "XI", 12: "XII", 13: "XIII", 14: "XIV", 15: "XV",
        16: "XVI", 17: "XVII", 18: "XVIII", 19: "XIX", 20: "XX",
    }
    roman_to_number = {value: key for key, value in roman.items()}

    for base in list(aliases):
        alternatives = [base]
        if re.search(r"\bdisease\b", base, flags=re.IGNORECASE):
            alternatives.append(
                re.sub(r"\bdisease\b", "syndrome", base, flags=re.IGNORECASE)
            )
        elif re.search(r"\bsyndrome\b", base, flags=re.IGNORECASE):
            alternatives.append(
                re.sub(r"\bsyndrome\b", "disease", base, flags=re.IGNORECASE)
            )

        for alternative in alternatives:
            add(alternative)
            type_match = re.search(
                r"\btype\s+([0-9]+|[IVXLCDM]+)([A-Z]?)\b",
                alternative,
                flags=re.IGNORECASE,
            )
            if type_match:
                number_text, suffix = type_match.groups()
                number = (
                    int(number_text)
                    if number_text.isdigit()
                    else roman_to_number.get(number_text.upper())
                )
                if number in roman:
                    for type_value in (str(number), roman[number]):
                        add(
                            alternative[:type_match.start(1)]
                            + type_value
                            + suffix.upper()
                            + alternative[type_match.end(2):]
                        )

            marker = re.match(
                r"^(.+?\b(?:disease|syndrome|disorder))\b",
                alternative,
                flags=re.IGNORECASE,
            )
            if marker:
                short_name = marker.group(1)
                add(short_name)
                if re.search(r"\bdisease\b", short_name, flags=re.IGNORECASE):
                    add(re.sub(
                        r"\bdisease\b", "syndrome", short_name,
                        flags=re.IGNORECASE,
                    ))
                elif re.search(r"\bsyndrome\b", short_name, flags=re.IGNORECASE):
                    add(re.sub(
                        r"\bsyndrome\b", "disease", short_name,
                        flags=re.IGNORECASE,
                    ))
    return aliases


class EuropePMCRequestError(RuntimeError):
    def __init__(self, message, partial=None):
        super().__init__(message)
        self.partial = partial


class InvalidEuropePMCResponse(ValueError):
    def __init__(self, message, partial=None):
        super().__init__(message)
        self.partial = partial


def validate_europe_pmc_response(payload, query):
    """Validate the first page before treating a response as complete or zero."""
    if not isinstance(payload, dict):
        raise InvalidEuropePMCResponse("Response must be a JSON object")
    raw_hits = payload.get("hitCount")
    valid_hits = not isinstance(raw_hits, bool) and (
        isinstance(raw_hits, int) or isinstance(raw_hits, str) and raw_hits.isdigit()
    )
    hits = int(raw_hits) if valid_hits else None
    if hits is not None and hits < 0:
        hits = None
    result_list = payload.get("resultList")
    raw_records = result_list.get("result") if isinstance(result_list, dict) else None
    if not isinstance(raw_records, list):
        raise InvalidEuropePMCResponse("Missing or invalid resultList.result")
    records = []
    for record in raw_records:
        if not isinstance(record, dict) or not isinstance(record.get("title"), str) or not record["title"].strip():
            continue
        pub_type_list = record.get("pubTypeList", {})
        if not isinstance(pub_type_list, dict):
            continue
        pub_types = pub_type_list.get("pubType", [])
        if isinstance(pub_types, str):
            pub_types = [pub_types]
        if not isinstance(pub_types, list) or not all(isinstance(value, str) for value in pub_types):
            continue
        records.append({
            "pmid": str(record.get("pmid", "")).strip(),
            "pmcid": str(record.get("pmcid", "")).strip(),
            "doi": str(record.get("doi", "")).strip(),
            "title": record["title"].strip(),
            "abstract": str(record.get("abstractText", "")).strip(),
            "year": str(record.get("pubYear") or record.get("firstPublicationDate") or "").strip(),
            "journal": str(record.get("journalTitle", "")).strip(),
            "cited_by": str(record.get("citedByCount", "")).strip(),
            "pub_types": pub_types,
        })
    partial = {"query": query, "hit_count": hits, "results": records}
    source_error = any(payload.get(field) for field in ("error", "errorCode", "errorMessage", "errCode", "errMsg"))
    if source_error or hits is None or len(records) != len(raw_records) or len(records) != min(hits, PAGE_SIZE):
        raise InvalidEuropePMCResponse("Invalid hitCount, record fields or incomplete first page", partial)
    return {"cache_version": CACHE_VERSION, "query": query, "hit_count": hits,
            "results": records, "response_validated": True}


def valid_cached_response(cached, query):
    if not isinstance(cached, dict) or cached.get("cache_version") != CACHE_VERSION or cached.get("query") != query:
        return False
    hits, records = cached.get("hit_count"), cached.get("results")
    if isinstance(hits, bool) or not isinstance(hits, int) or hits < 0 or not isinstance(records, list):
        return False
    if len(records) != min(hits, PAGE_SIZE):
        return False
    if not all(isinstance(record, dict) and isinstance(record.get("title"), str) and record["title"].strip()
               and isinstance(record.get("pub_types"), list) and all(isinstance(value, str) for value in record["pub_types"])
               for record in records):
        return False
    # Older positive caches preserve useful validated-shaped records. An older
    # zero could have been created from missing fields, so it needs revalidation.
    return hits > 0 or cached.get("response_validated") is True


def europe_pmc_request(query, cache_dir, session, counters):
    """Return one small cached Europe PMC response."""
    query = re.sub(r"\s+", " ", query).strip()
    cache_key = hashlib.sha256(
        json.dumps(
            [CACHE_VERSION, query, PAGE_SIZE, "core"],
            ensure_ascii=False,
        ).encode("utf-8")
    ).hexdigest()
    cache_path = cache_dir / f"{cache_key}.json"

    if cache_path.is_file():
        try:
            with open(cache_path, encoding="utf-8") as handle:
                cached = json.load(handle)
            if valid_cached_response(cached, query):
                counters["cache_hits"] += 1
                return cached
        except (OSError, json.JSONDecodeError):
            pass

    params = {
        "query": query,
        "format": "json",
        "resultType": "core",
        "pageSize": PAGE_SIZE,
    }
    last_error = None
    best_partial = None
    max_retries = getattr(session, "europe_pmc_max_retries", 3)
    for attempt in range(1, max_retries + 2):
        try:
            response = session.get(API_URL, params=params, timeout=60)
            response.raise_for_status()
            payload = response.json()
            result = validate_europe_pmc_response(payload, query)
            break
        except (requests.RequestException, ValueError) as error:
            last_error = error
            partial = getattr(error, "partial", None)
            if partial and partial["results"] and (best_partial is None or len(partial["results"]) >= len(best_partial["results"])):
                best_partial = partial
            if attempt == max_retries + 1:
                raise EuropePMCRequestError(f"Europe PMC request failed after {attempt} attempts: {query}", best_partial) from error
            time.sleep(2 ** (attempt - 1))
    else:  # defensive; the loop either succeeds or raises
        raise EuropePMCRequestError(f"Europe PMC request failed: {query}", best_partial) from last_error
    cache_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        suffix=".tmp",
        prefix=f".{cache_key}.",
        dir=cache_dir,
        delete=False,
    ) as temporary:
        json.dump(result, temporary, ensure_ascii=False, separators=(",", ":"))
        temporary_path = Path(temporary.name)
    try:
        os.replace(temporary_path, cache_path)
    finally:
        temporary_path.unlink(missing_ok=True)
    counters["api_requests"] += 1
    return result


def search_europe_pmc(
    gene,
    variant_identifiers,
    disease_terms,
    association_terms,
    cache_dir,
    session,
    counters,
):
    """Run field-restricted searches from specific evidence to rescue."""
    errors = []

    def request(query, *args):
        try:
            return europe_pmc_request(query, *args)
        except EuropePMCRequestError as error:
            errors.append(str(error))
            return error.partial or {"query": query, "hit_count": None, "results": []}

    def searchable_text(value, include_body=False):
        quoted = (
            '"'
            + str(value).replace("\\", "\\\\").replace('"', '\\"')
            + '"'
        )
        if include_body:
            return f"(TITLE_ABS:{quoted} OR BODY:{quoted})"
        return f"TITLE_ABS:{quoted}"

    def joined_terms(values, include_body=False):
        return " OR ".join(
            searchable_text(value, include_body) for value in values
        )

    def aliases_for(values):
        return list(dict.fromkeys(
            alias for value in values for alias in disease_aliases(value)
        ))

    def validated(response, required_groups=()):
        """Keep non-aggregate records with all terms in title or abstract."""
        if response is None:
            return None
        kept = []
        gene_pattern = re.compile(
            rf"(?<![A-Za-z0-9]){re.escape(gene)}(?![A-Za-z0-9])",
            flags=re.IGNORECASE,
        )
        aggregate_pattern = re.compile(
            r"^(?:abstracts?|conference abstracts?|meeting abstracts?|"
            r"e-?posters?|poster abstracts?|proceedings)(?:\b|\s|:)",
            flags=re.IGNORECASE,
        )
        for record in response["results"]:
            title = re.sub(r"<[^>]+>", " ", record.get("title", ""))
            abstract = re.sub(r"<[^>]+>", " ", record.get("abstract", ""))
            if aggregate_pattern.search(title.strip()):
                continue
            text = " " + re.sub(
                r"[^a-z0-9]+", " ", f"{title} {abstract}".lower()
            ).strip() + " "
            if not gene_pattern.search(f"{title} {abstract}"):
                continue
            if not all(
                any(
                    " "
                    + re.sub(r"[^a-z0-9]+", " ", term.lower()).strip()
                    + " " in text
                    for term in group
                )
                for group in required_groups
            ):
                continue
            kept.append(record)
        return {**response, "results": kept}

    exact_gene_query = searchable_text(gene, include_body=True)
    title_gene_query = searchable_text(gene)
    variant_query = joined_terms(variant_identifiers, include_body=True)
    disease_terms = aliases_for(disease_terms)
    disease_query = joined_terms(disease_terms)
    normalized_case_terms = {
        re.sub(r"\W+", " ", term.lower()).strip() for term in disease_terms
    }
    association_terms = [
        term for term in association_terms
        if re.sub(r"\W+", " ", term.lower()).strip() not in normalized_case_terms
    ]
    association_terms = aliases_for(association_terms)
    association_query = joined_terms(association_terms)
    layers = {
        "exact": None,
        "variant_context": None,
        "gene_context": None,
        "gene_association": None,
        "rescue": None,
    }

    if variant_query:
        layers["exact"] = request(
            f"{exact_gene_query} AND ({variant_query})",
            cache_dir,
            session,
            counters,
        )
    if variant_query and disease_query:
        layers["variant_context"] = request(
            f"{exact_gene_query} AND ({variant_query}) AND ({disease_query})",
            cache_dir,
            session,
            counters,
        )
    if disease_query:
        layers["gene_context"] = validated(
            request(
                f"{title_gene_query} AND ({disease_query})",
                cache_dir,
                session,
                counters,
            ),
            [disease_terms],
        )

    if association_query:
        counters["association_searches"] = counters.get("association_searches", 0) + 1
        layers["gene_association"] = validated(
            request(
                f"{title_gene_query} AND ({association_query})",
                cache_dir,
                session,
                counters,
            ),
            [association_terms],
        )

    retrieved = sum(
        len(value["results"]) for value in layers.values() if value is not None
    )
    if retrieved == 0:
        counters["rescue_searches"] += 1
        mutation_terms = ["mutation", "variant"]
        rare_disease_terms = ["rare disease", "genetic disorder"]
        layers["rescue"] = validated(
            request(
                f"{title_gene_query} AND ({joined_terms(mutation_terms)}) "
                f"AND ({joined_terms(rare_disease_terms)})",
                cache_dir,
                session,
                counters,
            ),
            [mutation_terms, rare_disease_terms],
        )
    if errors:
        layers["_errors"] = errors
    return layers

def summarize_publications(layers, study_patterns=None):
    """Deduplicate, classify and select representative publications."""
    if study_patterns is None:
        study_patterns = STUDY_PATTERNS

    layer_labels = [
        ("exact", "Exact variant publications retrieved"),
        ("variant_context", "Variant-context publications retrieved"),
        ("gene_context", "Gene-context publications retrieved"),
        ("gene_association", "Gene-association publications retrieved"),
    ]

    support_parts = []

    for key, label in layer_labels:
        value = layers.get(key)
        if value and value["results"]:
            support_parts.append(f"{label}: {len(value['results'])}")

    if layers.get("rescue") is not None and layers["rescue"]["results"]:
        support_parts.append(
            f"Rescue publications retrieved: "
            f"{len(layers['rescue']['results'])}"
        )

    unique = []
    seen_pmids, seen_dois, seen_titles = set(), set(), set()

    layer_order = (
        "exact",
        "variant_context",
        "gene_context",
        "gene_association",
        "rescue",
    )

    for layer in layer_order:
        value = layers.get(layer)
        if value is None:
            continue

        for record in value["results"]:
            pmid = str(record.get("pmid", "")).lower()
            doi = re.sub(
                r"^https?://(?:dx\.)?doi\.org/",
                "",
                str(record.get("doi", "")).lower(),
            )
            title = re.sub(
                r"\W+",
                " ",
                str(record.get("title", "")).lower(),
            ).strip()

            if (
                (pmid and pmid in seen_pmids)
                or (doi and doi in seen_dois)
                or (title and title in seen_titles)
            ):
                continue

            unique.append({
                "record": record,
                "layer": layer,
                "pmc_rank": len(unique),
                "study_types": set(),
            })

            if pmid:
                seen_pmids.add(pmid)
            if doi:
                seen_dois.add(doi)
            if title:
                seen_titles.add(title)

    reviewed = unique[:MAX_REVIEWED]
    if layers.get("_errors"):
        if not reviewed:
            # Exhausted searches cannot establish that the literature has zero
            # publications, even if some other layers returned valid zeros.
            return "", ""
        support_parts.insert(0, "Partial evidence: search incomplete")
    support_parts.append(
        f"Total publications reviewed: {len(reviewed)}"
    )

    study_counts = {name: 0 for name in study_patterns}

    for item in reviewed:
        record = item["record"]
        text = " ".join([
            str(record.get("title", "")),
            str(record.get("abstract", "")),
            " ".join(record.get("pub_types", [])),
        ])

        for name, pattern in study_patterns.items():
            if pattern.search(text):
                item["study_types"].add(name)
                study_counts[name] += 1

    detected = [
        f"{name}={count}"
        for name, count in study_counts.items()
        if count
    ]

    if reviewed:
        support_parts.append(
            "Study types detected: "
            + (
                ", ".join(detected)
                if detected
                else "not automatically classified"
            )
        )

    def publication_year(item):
        match = re.search(
            r"\b(?:19|20)\d{2}\b",
            str(item["record"].get("year", "")),
        )
        return int(match.group()) if match else -1

    selected = []
    selected_ranks = set()

    def add_best(predicate):
        candidates = [
            item
            for item in reviewed
            if item["pmc_rank"] not in selected_ranks
            and predicate(item)
        ]

        if not candidates:
            return False

        best = min(
            candidates,
            key=lambda item: (
                -publication_year(item),
                item["pmc_rank"],
            ),
        )

        selected.append(best)
        selected_ranks.add(best["pmc_rank"])
        return True

    # First select one variant-specific publication.
    if not add_best(lambda item: item["layer"] == "exact"):
        add_best(lambda item: item["layer"] == "variant_context")

    # Then add complementary study types.
    for study_type in ("human-cases", "functional", "animal"):
        if len(selected) >= MAX_TOP_PUBLICATIONS:
            break

        already_represented = any(
            study_type in item["study_types"]
            for item in selected
        )
        if not already_represented:
            add_best(
                lambda item, name=study_type:
                name in item["study_types"]
            )

    # Fill remaining positions by the original Europe PMC relevance order.
    for item in reviewed:
        if len(selected) >= MAX_TOP_PUBLICATIONS:
            break
        if item["pmc_rank"] not in selected_ranks:
            selected.append(item)
            selected_ranks.add(item["pmc_rank"])

    top = []

    for item in selected:
        record = item["record"]

        identifier = (
            f"PMID:{record['pmid']}" if record.get("pmid")
            else f"PMCID:{record['pmcid']}" if record.get("pmcid")
            else f"DOI:{record['doi']}" if record.get("doi")
            else "Europe PMC record"
        )

        citation = " | ".join(
            value
            for value in [
                identifier,
                record.get("year", ""),
                record.get("title", ""),
            ]
            if value
        )

        if record.get("doi") and not identifier.startswith("DOI:"):
            citation += f" | DOI:{record['doi']}"

        top.append(citation)

    return " | ".join(support_parts), " ;; ".join(top)

def main(argv=None):
    args = parse_args(argv)
    input_path = Path(args.input)
    output_path = Path(args.output)
    context_path = Path(args.case_context)
    cache_dir = Path(args.cache_dir)
    for path, label in ((input_path, "input"), (context_path, "case context")):
        if not path.is_file():
            raise FileNotFoundError(f"Stage 5 literature {label} not found: {path}")

    print("Stage 5 literature: loading input")
    input_df = pd.read_csv(input_path, dtype=str, keep_default_na=False)
    validate_table(input_df, "stage5_literature input", required=["MANE_Transcript_ID"], aligned=True)
    original_df = input_df.copy(deep=True)
    overlap = [column for column in OUTPUT_COLUMNS if column in input_df.columns]
    if overlap:
        raise ValueError("Literature columns already exist: " + ", ".join(overlap))

    def find_column(aliases, required=False):
        lower = {str(column).lower(): column for column in input_df.columns}
        for alias in aliases:
            if alias.lower() in lower:
                return lower[alias.lower()]
        if required:
            raise ValueError(f"Input is missing a required column: {aliases}")
        return None

    gene_col = find_column(GENE_COLUMN_ALIASES, required=True)
    aligned_col = ALIGNED_GENE_COLUMN if ALIGNED_GENE_COLUMN in input_df.columns else gene_col
    hgvsc_col = find_column(HGVSC_ALIASES)
    hgvsp_col = find_column(HGVSP_ALIASES)
    rsid_col = find_column(RSID_ALIASES)
    disease_associations_col = find_column(("Disease_associations",), required=True)
    literature_context = load_literature_context(context_path)
    disease_terms = literature_context["disease_terms"]

    print(f"  Rows: {len(input_df):,}")
    print(f"  Gene column: {gene_col}")
    print(f"  Transcript-aligned gene column: {aligned_col}")
    print(f"  Exact identifier columns: HGVSc={hgvsc_col}, HGVSp={hgvsp_col}, rsID={rsid_col}")
    print(f"  Gene-disease associations: {disease_associations_col}")
    print(
        f"  Disease context: {literature_context['context_source']} "
        f"({len(disease_terms)} terms)"
    )

    counters = {
        "api_requests": 0,
        "cache_hits": 0,
        "reused_units": 0,
        "association_searches": 0,
        "rescue_searches": 0,
        "gene_variant_units": 0,
    }
    source_errors = []
    failed_units = set()
    completed = {}
    support_values, top_values = [], []
    session = requests.Session()
    session.europe_pmc_max_retries = args.max_retries
    progress_every = max(0, args.progress_every)

    print("Stage 5 literature: searching Europe PMC")
    for row_number, (_, row) in enumerate(input_df.iterrows(), start=1):
        genes, identifiers = extract_gene_variants(
            row, gene_col, aligned_col, hgvsc_col, hgvsp_col, rsid_col
        )
        associated_diseases = extract_gene_diseases(
            row.get(disease_associations_col, ""), genes
        )
        row_support, row_top = {}, {}
        row_failed_genes = set()
        for gene in genes:
            counters["gene_variant_units"] += 1
            exact_terms = tuple(identifiers.get(gene, []))
            gene_diseases = tuple(associated_diseases.get(gene, []))
            search_key = (gene, exact_terms, tuple(disease_terms), gene_diseases)
            if search_key in completed:
                counters["reused_units"] += 1
                support, top = completed[search_key]
                if search_key in failed_units:
                    row_failed_genes.add(gene)
            else:
                try:
                    layers = search_europe_pmc(
                        gene,
                        exact_terms,
                        disease_terms,
                        gene_diseases,
                        cache_dir,
                        session,
                        counters,
                    )
                    support, top = summarize_publications(layers)
                    if layers.get("_errors"):
                        failed_units.add(search_key)
                        row_failed_genes.add(gene)
                        for error in layers["_errors"]:
                            warning = f"{gene}: Europe PMC review incomplete: {error}"
                            source_errors.append(warning)
                            print(f"Warning: {warning}")
                except Exception as error:  # source errors leave evidence blank
                    warning = f"{gene}: Europe PMC evidence unavailable: {error}"
                    source_errors.append(warning)
                    print(f"Warning: {warning}")
                    failed_units.add(search_key)
                    row_failed_genes.add(gene)
                    support, top = "", ""
                completed[search_key] = (support, top)
            row_support[gene] = support
            row_top[gene] = top

        if len(genes) > 1:
            support_values.append(" || ".join(f"{gene}: {row_support[gene]}" for gene in genes))
            top_values.append(" || ".join(
                f"{gene}: {row_top[gene] if gene in row_failed_genes else (row_top[gene] or 'No publications retrieved')}"
                for gene in genes
            ))
        elif genes:
            support_values.append(row_support[genes[0]])
            top_values.append(row_top[genes[0]])
        else:
            support_values.append("")
            top_values.append("")

        if progress_every and (row_number % progress_every == 0 or row_number == len(input_df)):
            print(
                f"  Rows {row_number:,}/{len(input_df):,} | "
                f"gene-variant units {counters['gene_variant_units']:,} | "
                f"API {counters['api_requests']:,} | cache {counters['cache_hits']:,} | "
                f"reused {counters['reused_units']:,} | "
                f"association {counters['association_searches']:,} | "
                f"rescue {counters['rescue_searches']:,}"
            )

    print("Stage 5 literature: assembling output")
    output_df = input_df.copy(deep=True)
    output_df["Literature_support"] = support_values
    output_df["Top_publications"] = top_values
    if len(output_df) != len(original_df):
        raise AssertionError("Stage 5 literature changed the number of rows")
    if not output_df[list(original_df.columns)].equals(original_df):
        raise AssertionError("Stage 5 literature changed existing values or row order")
    if list(output_df.columns) != list(original_df.columns) + OUTPUT_COLUMNS:
        raise AssertionError("Stage 5 literature produced unexpected columns")

    check_transition(original_df, output_df, "stage5_literature", aligned=True, warnings=source_errors)

    print("Stage 5 literature: writing output")
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
        temporary_path.unlink(missing_ok=True)

    print("\nStage 5 literature summary")
    print(f"Input rows:          {len(input_df):,}")
    print(f"Output rows:         {len(output_df):,}")
    print(f"Gene-variant units:  {counters['gene_variant_units']:,}")
    print(f"API requests:        {counters['api_requests']:,}")
    print(f"Cache hits:          {counters['cache_hits']:,}")
    print(f"Reused units:        {counters['reused_units']:,}")
    print(f"Association searches: {counters['association_searches']:,}")
    print(f"Rescue searches:     {counters['rescue_searches']:,}")


if __name__ == "__main__":
    main()
