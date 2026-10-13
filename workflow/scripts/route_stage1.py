import re
import pysam
from workflow_checks import VCFCounts, record_qc
import subprocess
from collections import Counter
from pathlib import Path

input_vcf = str(snakemake.input.vcf)
output_vcf = str(snakemake.output.vcf)
max_af = float(snakemake.params.max_af)

src = pysam.VariantFile(input_vcf)

# ---------------------------------------------------------
# Read CSQ structure from VEP header
# ---------------------------------------------------------

csq_desc = src.header.info["CSQ"].description
match = re.search(r"Format: (.+)", csq_desc)

if not match:
    raise RuntimeError("Could not parse CSQ field definition")

csq_fields = match.group(1).split("|")
csq_index = {name: i for i, name in enumerate(csq_fields)}

required = ["Consequence", "SYMBOL", "gnomADg_AF"]

for field in required:
    if field not in csq_index:
        raise RuntimeError(
            f"Missing CSQ field: {field}\n"
            f"Available fields: {csq_fields}"
        )

print("CSQ fields:", csq_fields)


# ---------------------------------------------------------
# Build output header
# ---------------------------------------------------------

header = src.header.copy()

header.info.add(
    "S1_SYMBOL",
    1,
    "String",
    "Gene symbol from Stage 1 minimal VEP annotation"
)

header.info.add(
    "S1_Consequence",
    1,
    "String",
    "Consequence from Stage 1 minimal VEP annotation"
)

header.info.add(
    "S1_gnomADg_AF",
    1,
    "Float",
    "gnomAD genome AF from VEP cache used for Stage 1 routing"
)

out = pysam.VariantFile(
    output_vcf,
    "wz",
    header=header
)


# ---------------------------------------------------------
# Route variants
# ---------------------------------------------------------

counts = VCFCounts(src.header.contigs)
annotated_by_contig = Counter()
n_total = 0
n_kept = 0
n_common = 0
n_missing = 0

for rec in src:
    counts.observe(rec.chrom, rec.pos, rec.ref, ",".join(rec.alts or []))
    annotated_by_contig[rec.chrom] += 1
    n_total += 1

    csq_values = rec.info.get("CSQ")

    symbol = None
    consequence = None
    af = None

    if csq_values:
        csq = csq_values[0]
        values = csq.split("|")
        if len(values) != len(csq_fields):
            raise ValueError("Malformed Stage 1 CSQ fields")

        symbol_raw = values[csq_index["SYMBOL"]]
        consequence_raw = values[csq_index["Consequence"]]
        af_raw = values[csq_index["gnomADg_AF"]]

        symbol = symbol_raw or None
        consequence = consequence_raw or None

        if af_raw:
            try:
                af = float(af_raw)
            except ValueError:
                af = None

    # Remove only known common variants
    if af is not None and af > max_af:
        n_common += 1
        continue

    if af is None:
        n_missing += 1

    # IMPORTANT:
    # Move the record from the original VEP header
    # to the new header containing S1_* definitions
    rec.translate(header)

    # We no longer need the full VEP CSQ payload
    if "CSQ" in rec.info:
        del rec.info["CSQ"]

    if symbol is not None:
        rec.info["S1_SYMBOL"] = symbol

    if consequence is not None:
        rec.info["S1_Consequence"] = consequence

    if af is not None:
        rec.info["S1_gnomADg_AF"] = af

    out.write(rec)
    n_kept += 1


src.close()
out.close()

source_stats = subprocess.check_output(
    ["bcftools", "index", "-s", str(snakemake.input.source)], text=True
)
source_by_contig = {}
for line in source_stats.splitlines():
    chrom, _length, number = line.split("\t")
    source_by_contig[chrom] = int(number)
source_count = sum(source_by_contig.values())

# Offline VEP may omit entire unplaced, unlocalized or decoy contigs that
# have no annotation features. Accept only these explicitly logged omissions;
# canonical chromosomes and every represented contig retain the count check.
omitted_by_contig = {}
for chrom in source_by_contig.keys() | annotated_by_contig.keys():
    expected = source_by_contig.get(chrom, 0)
    observed = annotated_by_contig.get(chrom, 0)
    if expected == observed:
        continue
    unsupported_name = (
        chrom.startswith("Un_")
        or chrom.endswith("_random")
        or chrom.endswith("_decoy")
    )
    if not unsupported_name or observed != 0 or expected <= 0:
        raise ValueError(
            f"Stage 1 annotation changed variant counts on {chrom}: "
            f"normalized={expected}, annotated={observed}"
        )
    omitted_by_contig[chrom] = expected

if omitted_by_contig:
    logged_skips = Counter()
    warning_path = Path(input_vcf + "_warnings.txt")
    if not warning_path.is_file():
        raise ValueError("Missing VEP evidence for unsupported-contig omissions")
    with warning_path.open() as warning_file:
        for line in warning_file:
            skipped = re.search(r"^WARNING: line \d+ skipped \((\S+) ", line)
            if (skipped and "not found in annotation sources or synonyms;" in line
                    and "does not overlap any features" in line):
                logged_skips[skipped.group(1)] += 1
    if dict(logged_skips) != omitted_by_contig:
        raise ValueError("VEP unsupported-contig warnings do not reconcile with omitted records")

if n_total != n_kept + n_common:
    raise ValueError("Stage 1 annotation/routing lost or duplicated variants")
omitted_count = sum(omitted_by_contig.values())
record_qc(
    "stage1_route", source_count, n_kept, source_count, n_kept, "filter",
    warnings=(["Offline VEP omitted documented unsupported noncanonical contigs"]
              if omitted_count else []),
    removed_variants=n_common + omitted_count,
    common_AF_removed=n_common,
    annotation_input_variants=n_total,
    unsupported_contig_variants_omitted=omitted_count,
    unsupported_contig_omissions=omitted_by_contig,
    missing_AF_kept=n_missing,
)

print()
print(f"Total variants:   {n_total:,}")
print(f"Unsupported contig variants omitted by VEP: {omitted_count:,}")
print(f"Removed AF >0.25: {n_common:,}")
print(f"Missing AF kept:  {n_missing:,}")
print(f"Stage 1 retained: {n_kept:,}")
