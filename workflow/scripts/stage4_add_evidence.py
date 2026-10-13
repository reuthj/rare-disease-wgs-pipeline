#!/usr/bin/env python3
"""Stage 4.2: add variant-level evidence without filtering Stage 4.1 rows.

The script keeps one row per normalized variant, preserves every Stage 4.1
column, checkpoints after each source, and resumes only when the checkpoint
matches the exact same input file and ordered variant keys.

Required environment: rare-wgs (pandas, pyarrow, pysam, pyBigWig).
SpliceAI is called through its installed command-line tool. ReMM v0.4 scores are read from the official remote Tabix
file; only its small index is kept in the technical Stage 4 directory.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
from workflow_checks import validate_table, check_transition, is_mitochondrial_variant


TX_SEP = " | "
NA = pd.NA

DEFAULT_INPUT = Path("results/intermediate/stage4/stage4_variants_routed.csv")
DEFAULT_OUTPUT = Path("results/tables/stage4_variant_evidence.csv")
DEFAULT_VEP_CACHE = Path("results/intermediate/stage3/vep_after_rarity.json.gz")
DEFAULT_STAGE4_DIR = Path("results/intermediate/stage4")
DEFAULT_GTF = Path("data/reference/annotation/gencode.v50.annotation.gtf.gz")
DEFAULT_CLINVAR_DIR = Path("data/resources/clinvar")
DEFAULT_ALPHAMISSENSE_DIR = Path("data/resources/alphamissense")

CHECKPOINT_NAME = "stage4_evidence_checkpoint.parquet"
STATE_NAME = "stage4_evidence_state.json"

ROUTE_COLUMNS = [
    "is_missense", "is_lof", "is_inframe_indel", "is_splice_relevant",
    "is_regulatory", "is_mitochondrial",
]
KEY_COLUMNS = ["Chr_Pos", "Variant_REF_ALT"]

EVIDENCE_COLUMNS = [
    "Known_variant",
    "ClinVar_conditions",
    "CADD",
    "phyloP100way",
    "phastCons100way",
    "Protein_context",
    "NMD_context",
    "Missense_evidence",
    "Splice_prediction",
    "Regulatory_prediction",
    "LoF_evidence",
    "Inframe_indel_evidence",
    "Mitochondrial_evidence",
]

REQUIRED_STAGE4_COLUMNS = {
    *KEY_COLUMNS,
    "Transcript_ID",
    "HGVSp",
    "Transcript_routes",
    "Evidence_routes",
    *ROUTE_COLUMNS,
}

PHYLOP_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/"
    "phyloP100way/hg38.phyloP100way.bw"
)
PHASTCONS_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/"
    "phastCons100way/hg38.phastCons100way.bw"
)
NCBI_DBSNP_BATCH_URL = (
    "https://api.ncbi.nlm.nih.gov/variation/v0/vcf/file/set_rsids"
)
NCBI_GRCH38_ASSEMBLY = "GCF_000001405.40"
REMM_URL = (
    "https://kircherlab.bihealth.org/download/ReMM/"
    "ReMM.v0.4.hg38.tsv.gz"
)
REMM_INDEX_MD5 = "c477b123abbe72c4aa66bf3abc18d9ab"


def banner(name: str) -> None:
    print("\n" + "=" * 60)
    print(name)
    print("=" * 60, flush=True)


def as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes"}


def has_value(value: Any) -> bool:
    if value is None or value is pd.NA:
        return False
    try:
        if pd.isna(value):
            return False
    except (TypeError, ValueError):
        pass
    return str(value).strip() not in {"", ".", "NA", "nan", "None"}


def clean_value(value: Any) -> str | None:
    return str(value).strip() if has_value(value) else None


def unique_strings(values: Iterable[Any]) -> list[str]:
    seen: set[str] = set()
    result: list[str] = []
    for value in values:
        text = clean_value(value)
        if text is not None and text not in seen:
            seen.add(text)
            result.append(text)
    return result


def normalize_chrom(value: Any) -> str:
    chrom = str(value).strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    return "MT" if chrom in {"M", "MT"} else chrom


def split_variant(row: pd.Series) -> tuple[str, int, str, str]:
    chrom, pos_text = str(row["Chr_Pos"]).split(":", 1)
    ref, alt = str(row["Variant_REF_ALT"]).split(">", 1)
    return normalize_chrom(chrom), int(pos_text), ref, alt


def variant_key(chrom: Any, pos: Any, ref: Any, alt: Any) -> str:
    return f"{normalize_chrom(chrom)}:{int(pos)}|{ref}>{alt}"


def row_key(row: pd.Series) -> str:
    return variant_key(*split_variant(row))


def parse_vep_key(result: dict[str, Any]) -> str | None:
    parts = str(result.get("input", "")).split()
    if len(parts) >= 5:
        return variant_key(parts[0], parts[1], parts[3], parts[4].split(",")[0])

    alleles = str(result.get("allele_string", "")).split("/")
    if len(alleles) == 2 and result.get("start") is not None:
        return variant_key(
            result.get("seq_region_name", ""),
            result["start"],
            alleles[0],
            alleles[1],
        )
    return None


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_keys(df: pd.DataFrame) -> str:
    digest = hashlib.sha256()
    for key in df.apply(row_key, axis=1):
        digest.update(key.encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def atomic_json_write(data: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    with temp.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, ensure_ascii=False)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temp, path)


def atomic_parquet_write(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    df.to_parquet(temp, index=False)
    os.replace(temp, path)


def atomic_csv_write(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    df.to_csv(temp, index=False, na_rep="NA")
    os.replace(temp, path)


def evidence_count(df: pd.DataFrame, column: str) -> int:
    return int(df[column].map(has_value).sum())


def route_count(df: pd.DataFrame, column: str) -> int:
    return int(df[column].map(as_bool).sum())


def transcript_map(value: Any) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for transcript in str(value or "").split(TX_SEP):
        transcript = transcript.strip()
        if not transcript or transcript == "-":
            continue
        mapping[transcript.split(".", 1)[0]] = transcript
    return mapping


def transcript_routes(value: Any) -> dict[str, set[str]]:
    mapping: dict[str, set[str]] = defaultdict(set)
    for item in str(value or "").split(TX_SEP):
        item = item.strip()
        if not item:
            continue
        separator = ":" if ":" in item else ("=" if "=" in item else None)
        if separator is None:
            continue
        transcript, routes = item.split(separator, 1)
        base = transcript.strip().split(".", 1)[0]
        for route in re.split(r"[+,;/|]", routes.upper()):
            route = route.strip()
            if route:
                mapping[base].add(route)
    return mapping


def format_tx_map(mapping: dict[str, list[str]]) -> Any:
    parts = []
    for transcript, values in mapping.items():
        clean = unique_strings(values)
        if clean:
            parts.append(f"{transcript}=" + ";".join(clean))
    return TX_SEP.join(parts) if parts else NA


def parse_tx_evidence(value: Any) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    if not has_value(value):
        return result
    for part in str(value).split(TX_SEP):
        if "=" not in part:
            continue
        transcript, evidence = part.split("=", 1)
        result[transcript] = [x for x in evidence.split(";") if x]
    return result


def add_tx_component(value: Any, transcript: str, component: str) -> Any:
    mapping = parse_tx_evidence(value)
    mapping.setdefault(transcript, [])
    if component not in mapping[transcript]:
        mapping[transcript].append(component)
    return format_tx_map(mapping)


def add_evidence_component(
    value: Any, component: str, *, prepend: bool = False
) -> str:
    components = unique_strings(str(value).split(TX_SEP)) if has_value(value) else []
    source = component.split("=", 1)[0]
    components = [
        item for item in components if item.split("=", 1)[0] != source
    ]
    if prepend:
        components.insert(0, component)
    else:
        components.append(component)
    return TX_SEP.join(components)


def format_number(value: Any) -> str | None:
    text = clean_value(value)
    if text is None:
        return None
    try:
        number = float(text)
        if math.isnan(number):
            return None
        return f"{number:.6g}"
    except ValueError:
        return text


def domains_text(value: Any) -> str | None:
    if not value:
        return None
    output = []
    entries = value if isinstance(value, list) else [value]
    for entry in entries:
        if isinstance(entry, dict):
            source = clean_value(entry.get("db") or entry.get("source"))
            name = clean_value(entry.get("name") or entry.get("id"))
            text = ":".join(x for x in [source, name] if x)
        else:
            text = clean_value(entry)
        if text:
            output.append(text)
    return ",".join(unique_strings(output)) or None


def protein_components(tc: dict[str, Any]) -> list[str]:
    components: list[str] = []
    amino = clean_value(tc.get("amino_acids"))
    if amino:
        components.append(f"AA={amino}")
    start = clean_value(tc.get("protein_start"))
    end = clean_value(tc.get("protein_end"))
    if start:
        position = start if not end or end == start else f"{start}-{end}"
        components.append(f"Protein_pos={position}")
    domains = domains_text(tc.get("domains"))
    if domains:
        components.append(f"Domains={domains}")
    return components


def consequence_supports(tc: dict[str, Any], route: str) -> bool:
    terms = {str(x).lower() for x in (tc.get("consequence_terms") or [])}
    if route == "MISSENSE":
        return "missense_variant" in terms
    if route == "LOF":
        return bool(
            terms
            & {
                "frameshift_variant",
                "stop_gained",
                "start_lost",
                "splice_acceptor_variant",
                "splice_donor_variant",
            }
        )
    if route == "INFRAME":
        return bool(
            terms & {"inframe_insertion", "inframe_deletion", "protein_altering_variant"}
        )
    return False


def transcript_supports(
    tc: dict[str, Any], route: str, route_map: dict[str, set[str]]
) -> bool:
    base = str(tc.get("transcript_id", "")).split(".", 1)[0]
    if route_map and base in route_map:
        return route in route_map[base]
    return consequence_supports(tc, route)


class Checkpoint:
    def __init__(self, directory: Path, input_path: Path, original: pd.DataFrame):
        self.checkpoint_path = directory / CHECKPOINT_NAME
        self.state_path = directory / STATE_NAME
        self.input_path = input_path
        self.input_hash = sha256_file(input_path)
        self.keys_hash = sha256_keys(original)
        self.original_columns = list(original.columns)

        checkpoint_exists = self.checkpoint_path.exists()
        state_exists = self.state_path.exists()
        if checkpoint_exists != state_exists:
            raise RuntimeError(
                "Checkpoint and state must either both exist or both be absent. "
                "Existing recovery files were not deleted automatically."
            )

        if checkpoint_exists:
            with self.state_path.open(encoding="utf-8") as handle:
                self.state = json.load(handle)
            self._verify_state(original)
            self.df = pd.read_parquet(self.checkpoint_path)
            self._verify_checkpoint(original)
            print(f"Resuming checkpoint: {self.checkpoint_path}")
            print("Completed steps: " + ", ".join(self.completed_steps))
        else:
            self.df = original.copy()
            for column in EVIDENCE_COLUMNS:
                if column not in self.df.columns:
                    self.df[column] = NA
            self.state = {
                "schema_version": 1,
                "input_path": str(input_path),
                "input_sha256": self.input_hash,
                "ordered_variant_keys_sha256": self.keys_hash,
                "input_rows": len(original),
                "original_columns": self.original_columns,
                "completed_steps": [],
                "step_summaries": {},
                "optional_sources_skipped": {},
            }

    @property
    def completed_steps(self) -> list[str]:
        return list(self.state.get("completed_steps", []))

    def _verify_state(self, original: pd.DataFrame) -> None:
        checks = {
            "input_sha256": self.input_hash,
            "ordered_variant_keys_sha256": self.keys_hash,
            "input_rows": len(original),
            "original_columns": self.original_columns,
        }
        for name, expected in checks.items():
            if self.state.get(name) != expected:
                raise RuntimeError(
                    f"Checkpoint does not match current input: {name}. "
                    "The existing checkpoint was left unchanged."
                )

    def _verify_checkpoint(self, original: pd.DataFrame) -> None:
        if len(self.df) != len(original):
            raise RuntimeError("Checkpoint row count differs from Stage 4.1 input")
        if sha256_keys(self.df) != self.keys_hash:
            raise RuntimeError("Checkpoint variant keys differ from Stage 4.1 input")
        missing = set(self.original_columns) - set(self.df.columns)
        if missing:
            raise RuntimeError(f"Checkpoint lost original columns: {sorted(missing)}")
        for column in ROUTE_COLUMNS:
            left = self.df[column].map(as_bool).tolist()
            right = original[column].map(as_bool).tolist()
            if left != right:
                raise RuntimeError(f"Checkpoint changed route column: {column}")

    def save(
        self,
        step: str,
        summary: dict[str, Any],
        optional_skip: str | None = None,
    ) -> None:
        if step not in self.state["completed_steps"]:
            self.state["completed_steps"].append(step)
        self.state["step_summaries"][step] = summary
        if optional_skip:
            self.state["optional_sources_skipped"][step] = optional_skip
        self.state["updated_at_epoch"] = time.time()

        # The dataframe is replaced first. If interruption occurs in the tiny
        # interval before state replacement, the prior state cannot claim that
        # data absent from the checkpoint was completed.
        atomic_parquet_write(self.df, self.checkpoint_path)
        atomic_json_write(self.state, self.state_path)
        print("Checkpoint saved", flush=True)

    def is_complete(self, step: str) -> bool:
        return step in self.completed_steps


def validate_initial_input(df: pd.DataFrame) -> None:
    validate_table(df, "Stage 4 evidence input", REQUIRED_STAGE4_COLUMNS | {"MANE_Transcript_ID"}, aligned=True)


def add_vep_evidence(df: pd.DataFrame, vep_cache: Path) -> dict[str, Any]:
    banner("VEP CACHE EVIDENCE")
    if not vep_cache.exists():
        raise FileNotFoundError(f"VEP cache not found: {vep_cache}")

    key_to_index = {row_key(row): index for index, row in df.iterrows()}
    found_keys: set[str] = set()

    cadd_by_index: dict[Any, dict[str, list[str]]] = defaultdict(dict)
    protein_by_index: dict[Any, dict[str, list[str]]] = defaultdict(dict)
    nmd_by_index: dict[Any, dict[str, list[str]]] = defaultdict(dict)
    missense_by_index: dict[Any, dict[str, list[str]]] = defaultdict(dict)
    lof_by_index: dict[Any, dict[str, list[str]]] = defaultdict(dict)
    inframe_by_index: dict[Any, dict[str, list[str]]] = defaultdict(dict)

    with gzip.open(vep_cache, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip().rstrip(",")
            if not line or line in {"[", "]"}:
                continue
            try:
                batch = json.loads(line)
            except json.JSONDecodeError as error:
                raise ValueError(f"Invalid VEP JSON at line {line_number}") from error

            if isinstance(batch, dict) and batch.get("complete") is False:
                raise RuntimeError(
                    f"VEP cache contains an incomplete batch at line {line_number}"
                )

            results = batch.get("results", []) if isinstance(batch, dict) else []
            for result in results:
                key = parse_vep_key(result)
                if key not in key_to_index:
                    continue
                index = key_to_index[key]
                found_keys.add(key)
                row = df.loc[index]
                tx_names = transcript_map(row["Transcript_ID"])
                routes = transcript_routes(row["Transcript_routes"])

                for tc in result.get("transcript_consequences", []) or []:
                    vep_tx = clean_value(tc.get("transcript_id"))
                    if not vep_tx:
                        continue
                    base = vep_tx.split(".", 1)[0]
                    if base not in tx_names:
                        continue
                    tx = tx_names[base]

                    cadd_parts = []
                    phred = format_number(tc.get("cadd_phred"))
                    raw = format_number(tc.get("cadd_raw"))
                    if phred:
                        cadd_parts.append(f"PHRED={phred}")
                    if raw:
                        cadd_parts.append(f"RAW={raw}")
                    if cadd_parts:
                        cadd_by_index[index][tx] = cadd_parts

                    protein = protein_components(tc)
                    coding_route = any(
                        transcript_supports(tc, route, routes)
                        for route in ["MISSENSE", "LOF", "INFRAME"]
                    )
                    if protein and coding_route:
                        protein_by_index[index][tx] = protein

                    if transcript_supports(tc, "MISSENSE", routes):
                        components = []
                        revel = format_number(tc.get("revel"))
                        sift_score = format_number(tc.get("sift_score"))
                        sift_prediction = clean_value(tc.get("sift_prediction"))
                        if revel:
                            components.append(f"REVEL={revel}")
                        if sift_score:
                            sift = f"SIFT={sift_score}"
                            if sift_prediction:
                                sift += f"({sift_prediction})"
                            components.append(sift)
                        if components:
                            missense_by_index[index][tx] = components

                    if transcript_supports(tc, "LOF", routes):
                        components = []
                        lof_value = clean_value(tc.get("lof"))
                        if lof_value:
                            components.append(f"LOFTEE={lof_value}")
                        for source, label in [
                            ("lof_filter", "FILTER"),
                            ("lof_flags", "FLAGS"),
                            ("lof_info", "INFO"),
                        ]:
                            value = clean_value(tc.get(source))
                            if value:
                                components.append(f"{label}={value}")
                        nmd_value = clean_value(tc.get("nmd"))
                        if nmd_value:
                            components.append(f"NMD={nmd_value}")
                            nmd_by_index[index][tx] = [nmd_value]
                        if components:
                            lof_by_index[index][tx] = components

                    if transcript_supports(tc, "INFRAME", routes) and protein:
                        inframe_by_index[index][tx] = protein

    for index in df.index:
        cadd = cadd_by_index.get(index, {})
        if cadd:
            unique_cadd = {tuple(values) for values in cadd.values()}
            if len(unique_cadd) == 1:
                df.at[index, "CADD"] = ";".join(next(iter(unique_cadd)))
            else:
                df.at[index, "CADD"] = format_tx_map(cadd)
        df.at[index, "Protein_context"] = format_tx_map(protein_by_index.get(index, {}))
        df.at[index, "NMD_context"] = format_tx_map(nmd_by_index.get(index, {}))
        df.at[index, "Missense_evidence"] = format_tx_map(
            missense_by_index.get(index, {})
        )
        df.at[index, "LoF_evidence"] = format_tx_map(lof_by_index.get(index, {}))
        df.at[index, "Inframe_indel_evidence"] = format_tx_map(
            inframe_by_index.get(index, {})
        )

    summary = {
        "eligible": len(df),
        "variants_found_in_cache": len(found_keys),
        "variants_not_found_in_cache": len(df) - len(found_keys),
        "CADD_annotated": evidence_count(df, "CADD"),
        "Missense_evidence_annotated": evidence_count(df, "Missense_evidence"),
        "LoF_evidence_annotated": evidence_count(df, "LoF_evidence"),
        "Protein_context_annotated": evidence_count(df, "Protein_context"),
    }
    for name, value in summary.items():
        print(f"{name}: {value:,}" if isinstance(value, int) else f"{name}: {value}")
    return summary


def post_with_retries(
    url: str,
    *,
    data: str,
    params: dict[str, str],
    attempts: int = 4,
) -> str:
    last_error: Exception | None = None
    query = urllib.parse.urlencode(params)
    request_url = f"{url}?{query}" if query else url
    for attempt in range(1, attempts + 1):
        try:
            request = urllib.request.Request(
                request_url,
                data=data.encode("utf-8"),
                headers={"Content-Type": "text/plain; charset=utf-8"},
                method="POST",
            )
            with urllib.request.urlopen(request, timeout=180) as response:
                return response.read().decode("utf-8")
        except (urllib.error.URLError, TimeoutError) as error:
            last_error = error
            if attempt < attempts:
                time.sleep(2 ** (attempt - 1))
    raise RuntimeError(f"NCBI dbSNP batch request failed: {last_error}")


def add_dbsnp(df: pd.DataFrame, *, _mitochondrial_batch: bool = False) -> dict[str, Any]:
    banner("dbSNP / KNOWN VARIANT")
    mitochondrial = df.apply(is_mitochondrial_variant, axis=1)
    if mitochondrial.any() and not _mitochondrial_batch:
        # A rejected MT lookup must not invalidate the nuclear batch.
        summaries = []
        warnings = []
        for is_mt in [False, True]:
            subset = df.loc[mitochondrial == is_mt].copy()
            if subset.empty:
                continue
            try:
                summaries.append(add_dbsnp(subset, _mitochondrial_batch=is_mt))
            except RuntimeError as error:
                if not is_mt:
                    raise
                warning = f"MT dbSNP evidence unavailable; existing identifiers retained: {error}"
                print(f"Warning: {warning}", flush=True)
                warnings.append(warning)
                continue
            df.loc[subset.index, "Known_variant"] = subset["Known_variant"]
        annotated = evidence_count(df, "Known_variant")
        return {"eligible": len(df), "annotated": annotated,
                "no_match": sum(summary["no_match"] for summary in summaries),
                "unavailable": int(mitochondrial.sum()) if warnings else 0,
                "warnings": warnings}
    lines = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    for _, row in df.iterrows():
        chrom, pos, ref, alt = split_variant(row)
        lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.")

    response_text = post_with_retries(
        NCBI_DBSNP_BATCH_URL,
        data="\n".join(lines) + "\n",
        params={"assembly": NCBI_GRCH38_ASSEMBLY},
        attempts=1 if _mitochondrial_batch else 4,
    )

    ids: dict[str, str] = {}
    response_errors: list[str] = []
    for line in response_text.splitlines():
        if not line:
            continue
        if line.startswith("# Error"):
            response_errors.append(line)
            continue
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 5:
            raise RuntimeError("NCBI returned malformed VCF output")
        if len(parts) >= 8 and "ERROR=" in parts[7]:
            response_errors.append(line)
            continue
        key = variant_key(parts[0], parts[1], parts[3], parts[4])
        rsid = parts[2].strip()
        if rsid not in {"", "."}:
            ids[key] = rsid

    if response_errors:
        raise RuntimeError(
            "NCBI dbSNP batch lookup returned variant-level errors; "
            f"first response lines: {response_errors[:4]}"
        )

    for index, row in df.iterrows():
        df.at[index, "Known_variant"] = ids.get(row_key(row), NA)

    annotated = evidence_count(df, "Known_variant")
    summary = {
        "eligible": len(df),
        "annotated": annotated,
        "no_match": len(df) - annotated,
    }
    print(f"Eligible variants: {len(df):,}")
    print(f"Annotated: {annotated:,}")
    print(f"No match: {len(df) - annotated:,}")
    return summary


def find_single_file(directory: Path, pattern: str, label: str) -> Path:
    matches = sorted(directory.rglob(pattern)) if directory.exists() else []
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected exactly one {label} matching {directory / pattern}; "
            f"found {len(matches)}"
        )
    return matches[0]


def flatten_info(value: Any) -> list[str]:
    values = value if isinstance(value, (tuple, list)) else [value]
    output = []
    for item in values:
        for part in str(item).split("|"):
            part = part.strip()
            if part and part not in {".", "not_provided"}:
                output.append(part)
    return output


def add_clinvar(df: pd.DataFrame, directory: Path) -> dict[str, Any]:
    banner("CLINVAR CONDITIONS")
    import pysam

    vcf_path = find_single_file(directory, "clinvar_*.vcf.gz", "ClinVar VCF")
    if not Path(str(vcf_path) + ".tbi").exists():
        raise FileNotFoundError(f"ClinVar Tabix index not found: {vcf_path}.tbi")

    vcf = pysam.VariantFile(str(vcf_path))
    contigs = set(vcf.header.contigs)
    annotated = 0

    try:
        for index, row in df.iterrows():
            chrom, pos, ref, alt = split_variant(row)
            candidates = [chrom, f"chr{chrom}"]
            if chrom == "MT":
                candidates.extend(["M", "chrM", "chrMT"])
            records = []
            for contig in candidates:
                if contig in contigs:
                    records.extend(vcf.fetch(contig, pos - 1, pos))

            evidence: list[str] = []
            for record in records:
                if record.pos != pos or record.ref != ref or alt not in (record.alts or ()):
                    continue
                conditions = flatten_info(record.info.get("CLNDN", ()))
                significance = flatten_info(record.info.get("CLNSIG", ()))
                if conditions and len(significance) == len(conditions):
                    evidence.extend(
                        f"{condition.replace('_', ' ')} = {sig.replace('_', ' ')}"
                        for condition, sig in zip(conditions, significance)
                    )
                elif conditions and significance:
                    sig_text = "/".join(unique_strings(significance)).replace("_", " ")
                    evidence.extend(
                        f"{condition.replace('_', ' ')} = {sig_text}"
                        for condition in conditions
                    )
                elif conditions:
                    evidence.extend(condition.replace("_", " ") for condition in conditions)

            clean = unique_strings(evidence)
            df.at[index, "ClinVar_conditions"] = TX_SEP.join(clean) if clean else NA
            annotated += bool(clean)
    finally:
        vcf.close()

    print(f"Eligible variants: {len(df):,}")
    print(f"Annotated: {annotated:,}")
    print(f"No match: {len(df) - annotated:,}")
    return {"eligible": len(df), "annotated": annotated, "no_match": len(df) - annotated}


def bigwig_value(
    bigwig: Any,
    available_chroms: set[str],
    chrom: str,
    pos: int,
) -> float | None:
    candidates = [f"chr{chrom}", chrom]
    if chrom == "MT":
        candidates = ["chrM", "chrMT", "MT", "M"]
    for candidate in candidates:
        if candidate not in available_chroms:
            continue
        values = bigwig.values(candidate, pos - 1, pos)
        if values and values[0] is not None and not math.isnan(values[0]):
            return float(values[0])
        return None
    return None


def bigwig_value_with_retries(
    pybigwig: Any,
    resource: dict[str, Any],
    chrom: str,
    pos: int,
    attempts: int = 5,
) -> float | None:
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            if resource["handle"] is None:
                handle = pybigwig.open(resource["url"])
                if handle is None:
                    raise RuntimeError(
                        f"Could not open remote BigWig: {resource['url']}"
                    )
                resource["handle"] = handle
                resource["chroms"] = set(handle.chroms())
            return bigwig_value(
                resource["handle"], resource["chroms"], chrom, pos
            )
        except Exception as error:
            last_error = error
            handle = resource.get("handle")
            if handle is not None:
                try:
                    handle.close()
                except Exception:
                    pass
            resource["handle"] = None
            resource["chroms"] = set()
            if attempt == attempts - 1:
                break
            delay = 2 ** (attempt + 1)
            print(
                f"{resource['column']}: remote read failed at "
                f"{chrom}:{pos}; retry {attempt + 2}/{attempts} "
                f"in {delay} seconds",
                flush=True,
            )
            time.sleep(delay)
    raise RuntimeError(
        f"Remote BigWig failed after {attempts} attempts at {chrom}:{pos}"
    ) from last_error


def add_conservation(df: pd.DataFrame) -> dict[str, Any]:
    banner("CONSERVATION")
    import pyBigWig

    resource_specs = [
        ("phyloP100way", PHYLOP_URL),
        ("phastCons100way", PHASTCONS_URL),
    ]
    resources = [
        {"column": column, "url": url, "handle": None, "chroms": set()}
        for column, url in resource_specs
    ]
    warnings = []
    try:
        for completed, (index, row) in enumerate(df.iterrows(), start=1):
            chrom, pos, _, _ = split_variant(row)
            for resource in resources:
                mitochondrial = is_mitochondrial_variant(row)
                try:
                    value = bigwig_value_with_retries(
                        pyBigWig, resource, chrom, pos,
                        attempts=1 if mitochondrial else 5,
                    )
                except RuntimeError as error:
                    if not mitochondrial:
                        raise
                    warning = f"{row_key(row)}: MT {resource['column']} unavailable: {error}"
                    print(f"Warning: {warning}", flush=True)
                    warnings.append(warning)
                    continue
                column = resource["column"]
                df.at[index, column] = (
                    format_number(value) if value is not None else NA
                )
            if completed % 100 == 0 or completed == len(df):
                print(
                    f"Conservation progress: {completed:,} / {len(df):,}",
                    flush=True,
                )
    except Exception as error:
        raise RuntimeError(
            "Remote UCSC BigWig lookup failed; no missing values were "
            "silently interpreted as successful no-match"
        ) from error
    finally:
        for resource in resources:
            handle = resource.get("handle")
            if handle is not None:
                handle.close()

    summary = {}
    for column, _ in resource_specs:
        annotated = evidence_count(df, column)
        summary[column] = {"annotated": annotated, "no_match": len(df) - annotated}
        print(f"{column}: {annotated:,} / {len(df):,}")
    if warnings:
        summary["warnings"] = warnings
    return summary


def tabix_header_columns(tabix: Any) -> list[str]:
    for line in reversed(list(tabix.header)):
        stripped = line.lstrip("#")
        if "\t" in stripped:
            return stripped.split("\t")
    raise ValueError("Could not find a tab-separated header in AlphaMissense file")


def resolve_column(columns: list[str], aliases: list[str], required: bool = True) -> str | None:
    lower = {column.lower().lstrip("#"): column for column in columns}
    for alias in aliases:
        if alias.lower() in lower:
            return lower[alias.lower()]
    if required:
        raise ValueError(f"None of the expected columns was found: {aliases}")
    return None


def add_alphamissense(df: pd.DataFrame, directory: Path) -> dict[str, Any]:
    banner("ALPHAMISSENSE")
    import pysam

    eligible_indices = [
        i for i in df.index
        if as_bool(df.at[i, "is_missense"]) and not is_mitochondrial_variant(df.loc[i])
    ]
    print(f"Eligible variants: {len(eligible_indices):,}")

    if not eligible_indices:
        return {
            "eligible": 0,
            "annotated": 0,
            "no_match": 0,
            "canonical": 0,
            "isoforms": 0,
        }

    source_specs = [
        (
            "canonical",
            "AlphaMissense_hg38.tsv.gz",
            "AlphaMissense canonical variant file",
        ),
        (
            "isoforms",
            "AlphaMissense_isoforms_hg38.tsv.gz",
            "AlphaMissense isoform file",
        ),
    ]

    sources: list[dict[str, Any]] = []

    try:
        # Open both indexed resources once.
        for label, filename, description in source_specs:
            path = find_single_file(directory, filename, description)
            index_path = Path(str(path) + ".tbi")

            if not index_path.exists():
                raise FileNotFoundError(
                    f"AlphaMissense indexed lookup is required, but the "
                    f"Tabix index is missing: {index_path}"
                )

            tabix = pysam.TabixFile(str(path))
            columns = tabix_header_columns(tabix)

            sources.append(
                {
                    "label": label,
                    "tabix": tabix,
                    "columns": columns,
                    "chrom_col": resolve_column(columns, ["CHROM", "chrom"]),
                    "pos_col": resolve_column(columns, ["POS", "position"]),
                    "ref_col": resolve_column(columns, ["REF", "ref"]),
                    "alt_col": resolve_column(columns, ["ALT", "alt"]),
                    "tx_col": resolve_column(
                        columns,
                        ["transcript_id", "transcript"],
                    ),
                    "score_col": resolve_column(
                        columns,
                        ["am_pathogenicity", "pathogenicity", "score"],
                    ),
                    "class_col": resolve_column(
                        columns,
                        [
                            "am_class",
                            "pathogenicity_class",
                            "classification",
                        ],
                        required=False,
                    ),
                }
            )

        annotated_variants = 0
        annotated_by_source = {
            "canonical": 0,
            "isoforms": 0,
        }

        for index in eligible_indices:
            row = df.loc[index]
            chrom, pos, ref, alt = split_variant(row)

            # Each transcript is matched independently.
            remaining_tx = dict(transcript_map(row["Transcript_ID"]))
            matched_variant = False

            if not remaining_tx:
                continue

            if str(chrom).startswith("chr"):
                contigs = [str(chrom)]
            else:
                contigs = [str(chrom), f"chr{chrom}"]

            # Canonical is checked first. Isoforms only fills transcripts
            # that were not matched in the canonical file.
            for source in sources:
                records: list[str] = []

                for contig in contigs:
                    try:
                        records.extend(
                            source["tabix"].fetch(contig, pos - 1, pos)
                        )
                    except (ValueError, OSError):
                        continue

                matched_source = False

                for line in records:
                    values = next(csv.reader([line], delimiter="\t"))
                    record = dict(zip(source["columns"], values))

                    if (
                        normalize_chrom(record[source["chrom_col"]]) != chrom
                        or int(record[source["pos_col"]]) != pos
                        or record[source["ref_col"]] != ref
                        or record[source["alt_col"]] != alt
                    ):
                        continue

                    record_tx = clean_value(
                        record.get(source["tx_col"])
                    )
                    if not record_tx:
                        continue

                    # ENST....10 and ENST....11 are treated as
                    # versions of the same stable transcript.
                    tx_base = record_tx.split(".", 1)[0]

                    if tx_base not in remaining_tx:
                        continue

                    score = format_number(
                        record.get(source["score_col"])
                    )
                    if not score:
                        continue

                    component = f"AlphaMissense={score}"

                    class_col = source["class_col"]
                    classification = (
                        clean_value(record.get(class_col))
                        if class_col
                        else None
                    )
                    if classification:
                        component += f"({classification})"

                    output_tx = remaining_tx.pop(tx_base)

                    df.at[index, "Missense_evidence"] = add_tx_component(
                        df.at[index, "Missense_evidence"],
                        output_tx,
                        component,
                    )

                    matched_source = True
                    matched_variant = True

                if matched_source:
                    annotated_by_source[source["label"]] += 1

                # All relevant transcripts have already been annotated.
                if not remaining_tx:
                    break

            annotated_variants += int(matched_variant)

    finally:
        for source in sources:
            source["tabix"].close()

    no_match = len(eligible_indices) - annotated_variants

    print(f"Annotated variants: {annotated_variants:,}")
    print(
        f"Canonical source: "
        f"{annotated_by_source['canonical']:,}"
    )
    print(
        f"Isoform fallback: "
        f"{annotated_by_source['isoforms']:,}"
    )
    print(f"No match: {no_match:,}")

    return {
        "eligible": len(eligible_indices),
        "annotated": annotated_variants,
        "no_match": no_match,
        "canonical": annotated_by_source["canonical"],
        "isoforms": annotated_by_source["isoforms"],
    }

def write_candidate_vcf(df: pd.DataFrame, indices: list[Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    records = [split_variant(df.loc[index]) for index in indices]
    with path.open("w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        for chrom in dict.fromkeys(record[0] for record in records):
            handle.write(f"##contig=<ID={chrom}>\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt in records:
            handle.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\n")


def find_reference_fasta() -> Path:
    preferred = Path(
        "data/reference/genome/GRCh38_no_alt_plus_hs38d1_no_chr.fasta"
    )
    if preferred.exists():
        return preferred
    matches = sorted(
        Path("data/reference").rglob("GRCh38_no_alt_plus_hs38d1_no_chr.fasta")
    )
    if len(matches) != 1:
        raise FileNotFoundError(
            "Could not resolve one GRCh38_no_alt_plus_hs38d1_no_chr.fasta "
            "under data/reference"
        )
    return matches[0]


def parse_gtf_attributes(text: str) -> dict[str, str]:
    result = {}
    for match in re.finditer(r'(\S+)\s+"([^"]*)"', text):
        result[match.group(1)] = match.group(2)
    return result


def build_spliceai_annotation(
    df: pd.DataFrame, indices: list[Any], gtf_path: Path, output: Path
) -> None:
    wanted: dict[str, str] = {}
    for index in indices:
        row = df.loc[index]
        for base, full in transcript_map(row["Transcript_ID"]).items():
            wanted[base] = full
        nearest = clean_value(row.get("Nearest_splice_transcript"))
        if nearest:
            wanted[nearest.split(".", 1)[0]] = nearest

    exons: dict[str, dict[str, Any]] = {}
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] != "exon":
                continue
            attrs = parse_gtf_attributes(parts[8])
            base = attrs.get("transcript_id", "").split(".", 1)[0]
            if base not in wanted:
                continue
            full = wanted[base]
            data = exons.setdefault(
                full,
                {"chrom": normalize_chrom(parts[0]), "strand": parts[6], "exons": []},
            )
            data["exons"].append((int(parts[3]) - 1, int(parts[4])))

    missing = sorted(set(wanted.values()) - set(exons))
    if missing:
        raise RuntimeError(
            f"GENCODE v50 is missing {len(missing)} requested splice transcripts; "
            f"first examples: {missing[:5]}"
        )

    temp = output.with_name(output.name + ".tmp")
    output.parent.mkdir(parents=True, exist_ok=True)
    with temp.open("w", encoding="utf-8") as handle:
        handle.write("#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END\n")
        for transcript, data in sorted(
            exons.items(), key=lambda item: (item[1]["chrom"], min(item[1]["exons"]))
        ):
            ordered = sorted(set(data["exons"]))
            starts = ",".join(str(start) for start, _ in ordered) + ","
            ends = ",".join(str(end) for _, end in ordered) + ","
            handle.write(
                f"{transcript}\t{data['chrom']}\t{data['strand']}\t"
                f"{ordered[0][0]}\t{ordered[-1][1]}\t{starts}\t{ends}\n"
            )
    os.replace(temp, output)


def parse_info_field(text: str) -> dict[str, str]:
    result = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def vcf_evidence(path: Path, info_name: str) -> dict[str, str]:
    evidence = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            value = parse_info_field(parts[7]).get(info_name)
            if value:
                evidence[variant_key(parts[0], parts[1], parts[3], parts[4])] = value
    return evidence


def run_command(
    command: list[str],
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
) -> None:
    print("Running: " + " ".join(str(x) for x in command), flush=True)
    subprocess.run(command, cwd=cwd, env=env, check=True)


def add_spliceai(
    df: pd.DataFrame,
    stage4_dir: Path,
    gtf_path: Path,
    conda_environment: str,
    conda_prefix: str | None = None,
) -> tuple[dict[str, Any], list[Path]]:
    banner("SPLICEAI")
    indices = [i for i in df.index if as_bool(df.at[i, "is_splice_relevant"]) and not is_mitochondrial_variant(df.loc[i])]
    print(f"Eligible variants: {len(indices):,}")
    if not indices:
        print("SPLICEAI: 0 variants — skipped")
        return {"eligible": 0, "annotated": 0, "no_match": 0}, []
    if not gtf_path.exists():
        raise FileNotFoundError(f"GENCODE v50 GTF not found: {gtf_path}")

    work = stage4_dir / "spliceai_work"
    work.mkdir(parents=True, exist_ok=True)
    input_vcf = work / "stage4_splice_candidates.vcf"
    output_vcf = work / "stage4_spliceai.vcf"
    annotation = work / "gencode.v50.stage4.spliceai.txt"
    write_candidate_vcf(df, indices, input_vcf)
    if not annotation.exists():
        build_spliceai_annotation(df, indices, gtf_path, annotation)

    reference = find_reference_fasta()
    if output_vcf.exists():
        print(f"Reusing completed technical SpliceAI output: {output_vcf}")
    else:
        executable, runtime_env = resolve_conda_runtime(
            conda_environment, "spliceai", environment_prefix=conda_prefix
        )
        temporary_output = work / f"stage4_spliceai.{os.getpid()}.tmp.vcf"
        run_command(
            [
                str(executable),
                "-I", str(input_vcf),
                "-O", str(temporary_output),
                "-R", str(reference),
                "-A", str(annotation),
                "-D", "50",
                "-M", "1",
            ],
            env=runtime_env,
        )
        if not temporary_output.exists():
            raise RuntimeError("SpliceAI completed without creating its output VCF")
        os.replace(temporary_output, output_vcf)

    raw = vcf_evidence(output_vcf, "SpliceAI")
    annotated = 0
    for index in indices:
        key = row_key(df.loc[index])
        entries = raw.get(key, "").split(",") if key in raw else []
        tx_names = transcript_map(df.at[index, "Transcript_ID"])
        mapping: dict[str, list[str]] = defaultdict(list)
        for entry in entries:
            fields = entry.split("|")
            if len(fields) < 10:
                continue
            transcript = fields[1]
            base = transcript.split(".", 1)[0]
            if base not in tx_names:
                continue
            scores = [float(x) for x in fields[2:6]]
            positions = fields[6:10]
            labels = ["AG", "AL", "DG", "DL"]
            maximum = max(scores)
            best = scores.index(maximum)
            mapping[tx_names[base]].append(
                f"SpliceAI={maximum:.6g};event={labels[best]};position={positions[best]}"
            )
        value = format_tx_map(mapping)
        df.at[index, "Splice_prediction"] = value
        annotated += has_value(value)

    print(f"Annotated: {annotated:,}")
    print(f"No match: {len(indices) - annotated:,}")
    return (
        {"eligible": len(indices), "annotated": annotated, "no_match": len(indices) - annotated},
        [work],
    )




def resolve_conda_runtime(
    environment_name: str, executable_name: str = "spliceai",
    environment_prefix: str | None = None,
) -> tuple[Path, dict[str, str]]:
    """Use an explicit prefix, the main runtime's sibling, or a unique named env."""
    if not environment_prefix and (not environment_name or Path(environment_name).name != environment_name):
        raise ValueError(
            "--spliceai-conda-env must be one Conda environment name, not a path"
        )
    sibling = Path(sys.prefix).resolve().parent / environment_name
    if environment_prefix:
        prefix = Path(environment_prefix).expanduser().resolve()
    elif sibling.is_dir():
        # Environments created together with the running rare-wgs environment
        # take precedence over older environments registered elsewhere.
        prefix = sibling.resolve()
    else:
        conda = os.environ.get("CONDA_EXE") or shutil.which("conda")
        if conda is None:
            raise RuntimeError("Conda executable was not found while resolving the predictor environment")
        result = subprocess.run(
            [conda, "env", "list", "--json"],
            check=True, capture_output=True, text=True,
        )
        try:
            prefixes = [Path(item) for item in json.loads(result.stdout).get("envs", [])]
        except (json.JSONDecodeError, AttributeError) as error:
            raise RuntimeError("Could not parse `conda env list --json` output") from error
        matches = sorted({path.resolve() for path in prefixes if path.name == environment_name})
        if len(matches) != 1:
            raise RuntimeError(
                f"Expected one Conda environment named {environment_name!r}; found {len(matches)}. "
                "Set stage4.spliceai_conda_prefix or --spliceai-conda-prefix to the intended environment."
            )
        prefix = matches[0]
    executable = prefix / "bin" / executable_name
    library_dir = prefix / "lib"
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise FileNotFoundError(
            f"{executable_name} executable not found: {executable}"
        )
    if not library_dir.is_dir():
        raise FileNotFoundError(f"Predictor Conda library directory not found: {library_dir}")
    print(f"SpliceAI environment: {prefix}", flush=True)

    runtime_env = os.environ.copy()
    runtime_env["CONDA_PREFIX"] = str(prefix)
    runtime_env["CONDA_DEFAULT_ENV"] = environment_name
    runtime_env["PATH"] = str(prefix / "bin") + os.pathsep + runtime_env.get("PATH", "")
    prior_library_path = runtime_env.get("LD_LIBRARY_PATH", "")
    runtime_env["LD_LIBRARY_PATH"] = str(library_dir) + (
        os.pathsep + prior_library_path if prior_library_path else ""
    )
    return executable, runtime_env




def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def ensure_remm_index(stage4_dir: Path) -> Path:
    """Keep only the small official index; the ~15 GB score file stays remote."""
    index_path = stage4_dir / "ReMM.v0.4.hg38.tsv.gz.tbi"
    if index_path.exists():
        if md5_file(index_path) != REMM_INDEX_MD5:
            raise RuntimeError(
                f"Existing ReMM index has an unexpected MD5: {index_path}. "
                "It was left unchanged."
            )
        return index_path

    stage4_dir.mkdir(parents=True, exist_ok=True)
    temporary = index_path.with_name(f"{index_path.name}.{os.getpid()}.tmp")
    request = urllib.request.Request(
        REMM_URL + ".tbi",
        headers={"User-Agent": "rare-wgs-stage4.2/1.0"},
    )
    try:
        with urllib.request.urlopen(request, timeout=180) as response:
            with temporary.open("wb") as output:
                shutil.copyfileobj(response, output)
                output.flush()
                os.fsync(output.fileno())
        observed_md5 = md5_file(temporary)
        if observed_md5 != REMM_INDEX_MD5:
            raise RuntimeError(
                "Downloaded ReMM v0.4 hg38 index failed its official MD5 check: "
                f"{observed_md5} != {REMM_INDEX_MD5}"
            )
        os.replace(temporary, index_path)
    finally:
        if temporary.exists():
            temporary.unlink()
    return index_path


def remm_positions(pos: int, ref: str, alt: str) -> list[int]:
    """Return the 1-based positions used by the official ReMM variant lookup."""
    if not ref or not alt or not re.fullmatch(r"[ACGTN]+", ref + alt, re.IGNORECASE):
        return []
    if len(ref) == len(alt):
        return list(range(pos, pos + len(ref)))
    if len(alt) > len(ref):
        return [pos, pos + 1]
    # A normalized VCF deletion includes the preceding anchor base in REF.
    return list(range(pos + 1, pos + len(ref)))


def contiguous_runs(positions: Iterable[int]) -> list[tuple[int, int]]:
    ordered = sorted(set(positions))
    if not ordered:
        return []
    runs: list[tuple[int, int]] = []
    start = previous = ordered[0]
    for current in ordered[1:]:
        if current != previous + 1:
            runs.append((start, previous))
            start = current
        previous = current
    runs.append((start, previous))
    return runs


def add_remm(df: pd.DataFrame, stage4_dir: Path) -> dict[str, Any]:
    banner("ReMM v0.4 (REMOTE TABIX)")
    import pysam

    indices = [i for i in df.index if as_bool(df.at[i, "is_regulatory"]) and not is_mitochondrial_variant(df.loc[i])]
    print(f"Eligible variants: {len(indices):,}")
    if not indices:
        print("ReMM: 0 variants — skipped")
        return {"eligible": 0, "annotated": 0, "no_match": 0, "version": "v0.4"}

    requested_by_row: dict[Any, tuple[str, list[int]]] = {}
    positions_by_chrom: dict[str, set[int]] = defaultdict(set)
    for index in indices:
        chrom, pos, ref, alt = split_variant(df.loc[index])
        positions = remm_positions(pos, ref, alt)
        requested_by_row[index] = (chrom, positions)
        positions_by_chrom[chrom].update(positions)

    index_path = ensure_remm_index(stage4_dir)
    tabix = None
    scores: dict[tuple[str, int], float] = {}
    try:
        tabix = pysam.TabixFile(REMM_URL, index=str(index_path))
        contig_by_chrom: dict[str, str] = {}
        for contig in tabix.contigs:
            contig_by_chrom.setdefault(normalize_chrom(contig), contig)

        for chrom, positions in positions_by_chrom.items():
            contig = contig_by_chrom.get(chrom)
            if contig is None:
                continue
            for start, end in contiguous_runs(positions):
                lines = None
                for attempt in range(5):
                    try:
                        lines = list(tabix.fetch(contig, start - 1, end))
                        break
                    except (OSError, ValueError) as error:
                        tabix.close()
                        tabix = None
                        if attempt == 4:
                            raise
                        delay = 2 ** (attempt + 1)
                        print(
                            f"ReMM remote read failed at {contig}:{start}-{end}; "
                            f"retrying in {delay} seconds",
                            flush=True,
                        )
                        time.sleep(delay)
                        tabix = pysam.TabixFile(REMM_URL, index=str(index_path))

                if lines is None:
                    raise RuntimeError(
                        f"ReMM returned no completed response for {contig}:{start}-{end}"
                    )

                for line in lines:
                    fields = line.split("\t")
                    if len(fields) < 3:
                        raise ValueError("Malformed row in ReMM Tabix response")
                    record_chrom = normalize_chrom(fields[0])
                    record_pos = int(fields[1])
                    if record_chrom != chrom or record_pos not in positions:
                        continue
                    score = float(fields[2])
                    if not 0.0 <= score <= 1.0:
                        raise ValueError(f"Invalid ReMM score at {fields[0]}:{fields[1]}")
                    key = (record_chrom, record_pos)
                    if key in scores and scores[key] != score:
                        raise ValueError(f"Conflicting ReMM scores at {fields[0]}:{fields[1]}")
                    scores[key] = score
    except Exception as error:
        raise RuntimeError(
            "Remote indexed ReMM lookup failed; missing values were not silently "
            "treated as successful no-match"
        ) from error
    finally:
        if tabix is not None:
            tabix.close()

    annotated = 0
    for index in indices:
        chrom, positions = requested_by_row[index]
        position_scores = [scores.get((chrom, position)) for position in positions]
        if not positions or any(score is None for score in position_scores):
            continue
        score = format_number(max(float(x) for x in position_scores if x is not None))
        if score is None:
            continue
        df.at[index, "Regulatory_prediction"] = add_evidence_component(
            df.at[index, "Regulatory_prediction"], f"ReMM={score}", prepend=True
        )
        annotated += 1

    print(f"Annotated: {annotated:,}")
    print(f"No match: {len(indices) - annotated:,}")
    return {
        "eligible": len(indices),
        "annotated": annotated,
        "no_match": len(indices) - annotated,
        "version": "v0.4",
    }














def mitochondrial_step(df: pd.DataFrame) -> dict[str, Any]:
    banner("MITOCHONDRIAL")
    eligible = int(df.apply(is_mitochondrial_variant, axis=1).sum())
    if eligible == 0:
        print("MITOCHONDRIAL: 0 variants — skipped")
        return {"eligible": 0, "annotated": 0, "no_match": 0, "initialized": False}
    warning = (
        "Mitochondrial-specific interpretation sources are not configured; "
        "unavailable evidence stays blank. Candidates continue with partial scoring."
    )
    print(f"Warning: {warning}", flush=True)
    return {"eligible": eligible, "annotated": 0, "no_match": 0,
            "initialized": False, "unavailable": eligible, "warnings": [warning]}


def cleanup_paths(paths: list[Path]) -> None:
    for path in paths:
        try:
            if path.is_dir():
                shutil.rmtree(path)
            elif path.exists():
                path.unlink()
        except OSError as error:
            print(f"Warning: could not remove technical temporary path {path}: {error}")


def source_step(
    result: tuple[dict[str, Any], list[Path]],
    optional_skip: str | None = None,
) -> tuple[dict[str, Any], str | None, list[Path]]:
    summary, cleanup = result
    return summary, optional_skip, cleanup


def validate_final(
    df: pd.DataFrame,
    original: pd.DataFrame,
) -> None:
    banner("FINAL VALIDATION")
    if len(df) != len(original):
        raise ValueError("Final row count changed")
    final_keys = df.apply(row_key, axis=1)
    original_keys = original.apply(row_key, axis=1)
    if final_keys.tolist() != original_keys.tolist():
        raise ValueError("A variant was added, removed, reordered, or changed")
    if final_keys.duplicated().any():
        raise ValueError("Final output contains duplicate normalized variant keys")
    missing_original = set(original.columns) - set(df.columns)
    if missing_original:
        raise ValueError(f"Original columns were lost: {sorted(missing_original)}")
    missing_evidence = set(EVIDENCE_COLUMNS) - set(df.columns)
    if missing_evidence:
        raise ValueError(f"Evidence columns are missing: {sorted(missing_evidence)}")
    for column in ROUTE_COLUMNS:
        before = original[column].map(as_bool).tolist()
        after = df[column].map(as_bool).tolist()
        if before != after:
            raise ValueError(f"Route column changed: {column}")

    print(f"Rows: {len(df):,}")
    print(f"Unique normalized variants: {final_keys.nunique():,}")
    print("\nEvidence coverage")
    for column in EVIDENCE_COLUMNS:
        print(f"{column:30s} {evidence_count(df, column):8,d}")
    print("\nFinal route counts")
    for column in ROUTE_COLUMNS:
        print(f"{column:30s} {route_count(df, column):8,d}")
    check_transition(original, df, "stage4_evidence", required=EVIDENCE_COLUMNS, aligned=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--vep-cache", type=Path, default=DEFAULT_VEP_CACHE)
    parser.add_argument("--stage4-dir", type=Path, default=DEFAULT_STAGE4_DIR)
    parser.add_argument("--gtf", type=Path, default=DEFAULT_GTF)
    parser.add_argument("--clinvar-dir", type=Path, default=DEFAULT_CLINVAR_DIR)
    parser.add_argument(
        "--alphamissense-dir", type=Path, default=DEFAULT_ALPHAMISSENSE_DIR
    )
    parser.add_argument(
        "--spliceai-conda-env",
        default="spliceai-env",
        help="Named Conda environment containing the SpliceAI executable",
    )
    parser.add_argument(
        "--spliceai-conda-prefix", default=None,
        help="Exact predictor environment directory (relative to project root or absolute); overrides name lookup",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    banner("STAGE 4.2 — VARIANT-LEVEL EVIDENCE")
    if not args.input.exists():
        raise FileNotFoundError(f"Exact Stage 4.1 input not found: {args.input}")

    original = pd.read_csv(args.input, dtype=str, keep_default_na=False)
    validate_initial_input(original)
    checkpoint = Checkpoint(args.stage4_dir, args.input, original)

    steps = [
        ("vep_cache", lambda: (add_vep_evidence(checkpoint.df, args.vep_cache), None, [])),
        ("dbsnp", lambda: (add_dbsnp(checkpoint.df), None, [])),
        ("clinvar", lambda: (add_clinvar(checkpoint.df, args.clinvar_dir), None, [])),
        ("conservation", lambda: (add_conservation(checkpoint.df), None, [])),
        (
            "alphamissense",
            lambda: (add_alphamissense(checkpoint.df, args.alphamissense_dir), None, []),
        ),
        (
            "spliceai",
            lambda: source_step(
                add_spliceai(
                    checkpoint.df,
                    args.stage4_dir,
                    args.gtf,
                    args.spliceai_conda_env,
                    args.spliceai_conda_prefix,
                )
            ),
        ),
        (
            "remm",
            lambda: (add_remm(checkpoint.df, args.stage4_dir), None, []),
        ),
        (
            "mitochondrial",
            lambda: (mitochondrial_step(checkpoint.df), None, []),
        ),
    ]

    for step, function in steps:
        if checkpoint.is_complete(step):
            print(f"{step}: completed in checkpoint — skipped", flush=True)
            continue
        result = function()
        summary, optional_skip, cleanup = result
        checkpoint.save(step, summary, optional_skip=optional_skip)
        cleanup_paths(cleanup)

    validate_final(checkpoint.df, original)
    atomic_csv_write(checkpoint.df, args.output)
    print(f"\nFinal output: {args.output}")


if __name__ == "__main__":
    main()
