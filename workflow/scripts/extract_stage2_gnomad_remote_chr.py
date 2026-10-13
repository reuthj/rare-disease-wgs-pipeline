#!/usr/bin/env python3

from datetime import datetime
from io import StringIO
from pathlib import Path
import os
import shlex
import signal
import subprocess
import time

import pandas as pd
from workflow_checks import record_qc


# ============================================================
# Snakemake inputs
# ============================================================

CHROM = str(snakemake.wildcards.chrom)
KEYS_BCF = str(Path(snakemake.input.variant_keys).resolve())
GNOMAD_URL = str(snakemake.params.url)
OUT = str(snakemake.output[0])
REMOTE_RETRIES = int(snakemake.params.retries)
RETRY_DELAY_SECONDS = float(snakemake.params.retry_delay_seconds)
REMOTE_TIMEOUT_SECONDS = float(getattr(snakemake.params, "remote_timeout_seconds", 7200))
INDEX_CACHE_DIR = Path(getattr(snakemake.params, "index_cache_dir",
                             "data/resources/gnomad/joint/v4.1/indexes")).resolve()


def run_remote(command, failure_context, timeout_seconds=None):
    """Bound remote work and reap the entire query pipeline before retrying."""
    timeout_seconds = REMOTE_TIMEOUT_SECONDS if timeout_seconds is None else timeout_seconds
    INDEX_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               text=True, start_new_session=True, cwd=INDEX_CACHE_DIR)
    try:
        stdout, stderr = process.communicate(timeout=timeout_seconds)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        try:
            process.communicate(timeout=10)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.communicate()
        raise RuntimeError(f"{failure_context}: timed out after {timeout_seconds:g} seconds") from None
    return subprocess.CompletedProcess(command, process.returncode, stdout, stderr)

# GRCh38 chrX pseudoautosomal regions
PAR1 = (10001, 2781479)
PAR2 = (155701383, 156030895)

KEY = ["_chrom", "_pos", "_ref", "_alt"]

# Only raw Joint fields needed to build the compact final Stage-2 fields.
RAW_FIELDS = [
    "AF_joint",
    "grpmax_joint",
    "AF_grpmax_joint",
    "AC_joint",
    "AN_joint",
    "nhomalt_joint",
    "AC_joint_XX",
    "nhomalt_joint_XX",
    "AC_joint_XY",
]

FINAL_FIELDS = [
    "gnomAD_AF",
    "gnomAD_grpmax_AF",
    "gnomAD_grpmax_population",
    "gnomAD_AC",
    "gnomAD_AN",
    "gnomAD_het",
    "gnomAD_hom",
    "gnomAD_hemi",
]


# ============================================================
# Helpers
# ============================================================

def log(message):
    print(
        f"[{datetime.now():%Y-%m-%d %H:%M:%S}] "
        f"[chr{CHROM}] {message}",
        flush=True,
    )


def clean_chrom(value):
    return str(value).removeprefix("chr")


def in_x_par(pos):
    return (
        PAR1[0] <= pos <= PAR1[1]
        or PAR2[0] <= pos <= PAR2[1]
    )


def expected_variant_count(contig):
    """
    Count project variants for this chromosome in the small indexed key BCF.
    """
    command = (
        f"bcftools view -r {shlex.quote(contig)} -H {shlex.quote(KEYS_BCF)} "
        "| wc -l"
    )

    result = subprocess.run(
        ["bash", "-o", "pipefail", "-c", command],
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to count target variants for {contig}:\n"
            f"{result.stderr}"
        )

    return int(result.stdout.strip())


def check_remote_header():
    """
    Verify that this gnomAD Joint file exposes the fields required by
    the compact Stage-2 annotation.
    """
    result = run_remote(
        ["bcftools", "view", "-h", GNOMAD_URL],
        f"Could not read gnomAD header for chr{CHROM}",
        timeout_seconds=min(120, REMOTE_TIMEOUT_SECONDS),
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"Could not read gnomAD header for chr{CHROM}:\n"
            f"{result.stderr}"
        )

    header = result.stdout
    missing = [
        tag
        for tag in RAW_FIELDS
        if f"##INFO=<ID={tag}," not in header
    ]

    if missing:
        raise ValueError(
            f"chr{CHROM}: gnomAD Joint file is missing required INFO fields: "
            f"{missing}"
        )


# ============================================================
# Remote gnomAD query
# ============================================================

def query_gnomad():
    """
    Annotate this chromosome's project variants directly from the remote
    gnomAD Joint v4.1 VCF.

    Core strategy kept from the previously successful script:
      indexed local BCF -> bcftools annotate -a remote_url -> bcftools query
    """
    contig = f"chr{CHROM}"
    expected_n = expected_variant_count(contig)

    if expected_n == 0:
        raise ValueError(
            f"chr{CHROM}: no project variants found in {KEYS_BCF}"
        )

    log(
        f"{expected_n:,} project variants; querying remote gnomAD..."
    )

    source_columns = ",".join(
        f"INFO/{tag}" for tag in RAW_FIELDS
    )

    fmt = (
        "%CHROM\\t"
        "%POS\\t"
        "%REF\\t"
        "%ALT\\t"
        + "\\t".join(f"%INFO/{tag}" for tag in RAW_FIELDS)
        + "\\n"
    )

    names = [
        "CHROM",
        "POS",
        "REF",
        "ALT",
        *RAW_FIELDS,
    ]

    command = f"""
    set -o pipefail

    bcftools annotate \
        -r {shlex.quote(contig)} \
        -a {shlex.quote(GNOMAD_URL)} \
        -c {shlex.quote(source_columns)} \
        {shlex.quote(KEYS_BCF)} \
        -Ou \
    | bcftools query \
        -f '{fmt}'
    """

    result = run_remote(
        ["bash", "-c", command],
        f"chr{CHROM} remote gnomAD query failed",
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"chr{CHROM} remote gnomAD query failed:\n"
            f"{result.stderr}"
        )

    try:
        df = pd.read_csv(
            StringIO(result.stdout),
            sep="\t",
            header=None,
            names=names,
            na_values=".",
            keep_default_na=True,
        )
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns=names)

    if len(df) != expected_n:
        raise ValueError(
            f"chr{CHROM}: expected {expected_n:,} project variants "
            f"but received {len(df):,} rows after annotation"
        )

    return df


# ============================================================
# Convert raw Joint fields -> compact Stage-2 fields
# ============================================================

def build_checkpoint(df):
    if df.empty:
        raise ValueError(f"chr{CHROM}: gnomAD query returned no variants")

    numeric_columns = [
        "POS",
        "AF_joint",
        "AF_grpmax_joint",
        "AC_joint",
        "AN_joint",
        "nhomalt_joint",
        "AC_joint_XX",
        "nhomalt_joint_XX",
        "AC_joint_XY",
    ]

    for column in numeric_columns:
        df[column] = pd.to_numeric(
            df[column],
            errors="coerce",
        )

    df["POS"] = df["POS"].astype("int64")

    for column in [
        "AC_joint",
        "AN_joint",
        "nhomalt_joint",
        "AC_joint_XX",
        "nhomalt_joint_XX",
        "AC_joint_XY",
    ]:
        df[column] = df[column].astype("Int64")

    df["_chrom"] = df["CHROM"].map(clean_chrom)
    df["_pos"] = df["POS"]
    df["_ref"] = df["REF"].astype(str)
    df["_alt"] = df["ALT"].astype(str)

    wrong_chrom = df["_chrom"] != CHROM
    if wrong_chrom.any():
        examples = df.loc[
            wrong_chrom,
            ["CHROM", "POS", "REF", "ALT"],
        ].head(5)

        raise ValueError(
            f"chr{CHROM}: output contains records from another chromosome:\n"
            + examples.to_string(index=False)
        )

    if df.duplicated(KEY).any():
        raise ValueError(
            f"chr{CHROM}: duplicate project variant keys after gnomAD annotation"
        )

    found = df["AC_joint"].notna()

    # Canonical compact Stage-2 aliases.
    df["gnomAD_AF"] = df["AF_joint"]
    df["gnomAD_grpmax_AF"] = df["AF_grpmax_joint"]
    df["gnomAD_grpmax_population"] = df["grpmax_joint"].astype("string")
    df["gnomAD_AC"] = df["AC_joint"].copy()
    df["gnomAD_AN"] = df["AN_joint"].copy()

    df["gnomAD_het"] = pd.Series(
        pd.NA,
        index=df.index,
        dtype="Int64",
    )
    df["gnomAD_hom"] = pd.Series(
        pd.NA,
        index=df.index,
        dtype="Int64",
    )
    df["gnomAD_hemi"] = pd.Series(
        pd.NA,
        index=df.index,
        dtype="Int64",
    )

    # --------------------------------------------------------
    # Autosomes
    # --------------------------------------------------------

    if CHROM not in {"X", "Y"}:
        valid = (
            found
            & df["AC_joint"].notna()
            & df["nhomalt_joint"].notna()
        )

        df.loc[valid, "gnomAD_hom"] = (
            df.loc[valid, "nhomalt_joint"]
        )

        df.loc[valid, "gnomAD_het"] = (
            df.loc[valid, "AC_joint"]
            - 2 * df.loc[valid, "nhomalt_joint"]
        )

        df.loc[valid, "gnomAD_hemi"] = 0

    # --------------------------------------------------------
    # chrX
    # --------------------------------------------------------

    elif CHROM == "X":
        par_mask = df["_pos"].map(
            lambda pos: in_x_par(int(pos))
        )

        # PAR behaves diploid.
        par = (
            found
            & par_mask
            & df["AC_joint"].notna()
            & df["nhomalt_joint"].notna()
        )

        df.loc[par, "gnomAD_hom"] = (
            df.loc[par, "nhomalt_joint"]
        )

        df.loc[par, "gnomAD_het"] = (
            df.loc[par, "AC_joint"]
            - 2 * df.loc[par, "nhomalt_joint"]
        )

        df.loc[par, "gnomAD_hemi"] = 0

        # Outside PAR: XX contributes diploid carriers, XY hemizygous alleles.
        nonpar = (
            found
            & ~par_mask
            & df["AC_joint_XX"].notna()
            & df["nhomalt_joint_XX"].notna()
            & df["AC_joint_XY"].notna()
        )

        df.loc[nonpar, "gnomAD_hom"] = (
            df.loc[nonpar, "nhomalt_joint_XX"]
        )

        df.loc[nonpar, "gnomAD_het"] = (
            df.loc[nonpar, "AC_joint_XX"]
            - 2 * df.loc[nonpar, "nhomalt_joint_XX"]
        )

        df.loc[nonpar, "gnomAD_hemi"] = (
            df.loc[nonpar, "AC_joint_XY"]
        )

    # --------------------------------------------------------
    # chrY
    # --------------------------------------------------------

    else:
        # Y is not treated as a diploid autosome.
        # het/hom are not applicable; hemi is the XY alternate count.
        valid = found & df["AC_joint_XY"].notna()

        df.loc[valid, "gnomAD_hemi"] = (
            df.loc[valid, "AC_joint_XY"]
        )

    # --------------------------------------------------------
    # Variants absent from gnomAD Joint
    # --------------------------------------------------------

    absent = ~found

    # No Joint row -> zero observed alternate alleles.
    # AN/AF/grpmax remain missing because no site row provides them.
    df.loc[absent, "gnomAD_AC"] = 0

    if CHROM == "Y":
        df.loc[absent, "gnomAD_hemi"] = 0
    else:
        df.loc[
            absent,
            ["gnomAD_het", "gnomAD_hom", "gnomAD_hemi"],
        ] = 0

    # Preserve nullable integer types.
    for column in [
        "gnomAD_AC",
        "gnomAD_AN",
        "gnomAD_het",
        "gnomAD_hom",
        "gnomAD_hemi",
    ]:
        df[column] = pd.to_numeric(
            df[column],
            errors="coerce",
        ).astype("Int64")

    # --------------------------------------------------------
    # Sanity checks
    # --------------------------------------------------------

    for column in ["gnomAD_AF", "gnomAD_grpmax_AF"]:
        bad = (
            df[column].notna()
            & (
                (df[column] < 0)
                | (df[column] > 1)
            )
        )

        if bad.any():
            raise ValueError(
                f"chr{CHROM}: {column} contains values outside [0, 1]"
            )

    for column in [
        "gnomAD_AC",
        "gnomAD_AN",
        "gnomAD_het",
        "gnomAD_hom",
        "gnomAD_hemi",
    ]:
        bad = df[column].notna() & (df[column] < 0)

        if bad.any():
            raise ValueError(
                f"chr{CHROM}: negative values found in {column}"
            )

    return df[KEY + FINAL_FIELDS].copy()


# ============================================================
# Main
# ============================================================

def query_with_retries():
    """Retry transport failures; scientific/schema validation remains fatal."""
    attempts = []
    for attempt in range(1, REMOTE_RETRIES + 2):
        log(f"remote attempt {attempt}/{REMOTE_RETRIES + 1}")
        started = time.monotonic()
        try:
            check_remote_header()
            raw = query_gnomad()
        except RuntimeError as error:
            message = str(error)
            if ("Could not read gnomAD header" not in message
                    and "remote gnomAD query failed" not in message):
                raise
            attempts.append(dict(attempt=attempt, status="FAILED",
                                 seconds=round(time.monotonic() - started, 3),
                                 error=message))
            log(f"remote attempt {attempt} failed: {message}")
            if attempt == REMOTE_RETRIES + 1:
                raise RuntimeError(
                    f"chr{CHROM}: remote gnomAD failed after {attempt} attempts "
                    f"(initial attempt plus {REMOTE_RETRIES} retries):\n{message}"
                ) from error
            delay = RETRY_DELAY_SECONDS * 2 ** (attempt - 1)
            log(f"retry {attempt}/{REMOTE_RETRIES} in {delay:g} seconds")
            time.sleep(delay)
        else:
            attempts.append(dict(attempt=attempt, status="PASS",
                                 seconds=round(time.monotonic() - started, 3)))
            return raw, attempts


def main():
    start = time.time()

    log("starting remote gnomAD annotation")

    expected_n = expected_variant_count(f"chr{CHROM}")
    attempts = []
    if expected_n == 0:
        out = pd.DataFrame(columns=KEY + FINAL_FIELDS)
    else:
        raw, attempts = query_with_retries()
        out = build_checkpoint(raw)
    if len(out) != expected_n:
        raise ValueError("Chromosome checkpoint changed case target counts")

    os.makedirs(
        os.path.dirname(OUT),
        exist_ok=True,
    )

    tmp = OUT + ".tmp"
    out.to_parquet(
        tmp,
        index=False,
    )
    os.replace(tmp, OUT)

    record_qc(f"stage2_gnomad_chr{CHROM}", expected_n, len(out), expected_n, len(out), "enrichment",
              configured_retries=REMOTE_RETRIES, remote_attempts=attempts,
              retries_used=max(0, len(attempts) - 1))
    elapsed = time.time() - start
    observed = int((out["gnomAD_AC"].fillna(0) > 0).sum())

    log(
        f"DONE | targets={len(out):,} | "
        f"observed_AC>0={observed:,} | "
        f"time={elapsed / 60:.1f} min | "
        f"saved={OUT}"
    )


if __name__ == "__main__":
    main()
