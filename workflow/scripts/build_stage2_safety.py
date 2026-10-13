import os
import re
import subprocess
from workflow_checks import VCFCounts, record_qc
import pysam

import pyarrow as pa
import pyarrow.parquet as pq


# Must match --fields in stage2_vep_context
CSQ_FIELDS = [
    "Allele",
    "Consequence",
    "SYMBOL",
    "Gene",
    "Feature_type",
    "Feature",
    "BIOTYPE",
    "EXON",
    "INTRON",
    "MANE_SELECT",
    "PICK",
    "ENCODE_cCRE",
]


def clean(value):
    if value in (None, "", "."):
        return None
    return value


def parse_csq(csq):
    csq = clean(csq)

    if csq is None:
        return []

    records = []

    for item in csq.split(","):
        values = item.split("|")

        if len(values) != len(CSQ_FIELDS):
            raise ValueError(f"Malformed Stage 2 CSQ: expected {len(CSQ_FIELDS)} fields, found {len(values)}")

        records.append(dict(zip(CSQ_FIELDS, values)))

    return records


def unique(values):
    result = []
    seen = set()

    for value in values:
        value = clean(value)

        if value is None or value in seen:
            continue

        seen.add(value)
        result.append(value)

    return result


def genomic_region(record):
    if not record:
        return None

    exon = clean(record.get("EXON"))
    intron = clean(record.get("INTRON"))

    if exon:
        return f"Exon {exon.split('/')[0]}"

    if intron:
        return f"Intron {intron.split('/')[0]}"

    consequences = set(
        clean(record.get("Consequence") or "").split("&")
    )

    if "5_prime_UTR_variant" in consequences:
        return "5'UTR"

    if "3_prime_UTR_variant" in consequences:
        return "3'UTR"

    if "upstream_gene_variant" in consequences:
        return "Upstream"

    if "downstream_gene_variant" in consequences:
        return "Downstream"

    if "intergenic_variant" in consequences:
        return "Intergenic"

    return None


def regulatory_context(records):
    ensembl = []
    encode = []

    for record in records:

        if record.get("Feature_type") == "RegulatoryFeature":
            biotype = clean(record.get("BIOTYPE"))
            feature = clean(record.get("Feature"))

            value = "Ensembl"

            if biotype:
                value += f": {biotype}"

            if feature:
                value += f" | {feature}"

            ensembl.append(value)

        ccre = clean(record.get("ENCODE_cCRE"))

        if ccre:
            encode.append(f"ENCODE: {ccre}")

    values = unique(ensembl + encode)

    if not values:
        return None

    return " ; ".join(values)


def selected_transcripts(records):
    transcripts = [
        r for r in records
        if r.get("Feature_type") == "Transcript"
    ]

    picked = [
        r for r in transcripts
        if r.get("PICK") == "1"
    ]

    if picked:
        return picked

    # Defensive fallback: one transcript per gene
    by_gene = {}

    for record in transcripts:
        key = (
            clean(record.get("Gene"))
            or clean(record.get("SYMBOL"))
            or "__unknown__"
        )

        if key not in by_gene:
            by_gene[key] = record

    return list(by_gene.values())


schema = pa.schema([
    ("Chr_Pos", pa.string()),
    ("Variant_REF_ALT", pa.string()),
    ("Gene_symbol", pa.string()),
    ("Gene_ID", pa.string()),
    ("Transcript_ID", pa.string()),
    ("MANE_Select", pa.string()),
    ("Transcript_biotype", pa.string()),
    ("Consequence", pa.string()),
    ("Genomic_region", pa.string()),
    ("ClinVar_ID", pa.string()),
    ("ClinVar_significance", pa.string()),
    ("Regulatory_context", pa.string()),
])


INPUT_VCF = str(snakemake.input.vcf)
OUTPUT_PARQUET = str(snakemake.output.parquet)

os.makedirs(
    os.path.dirname(OUTPUT_PARQUET),
    exist_ok=True
)


with pysam.VariantFile(INPUT_VCF) as header_file:
    fields_match = re.search(r"Format: (.+)", header_file.header.info["CSQ"].description)
    if not fields_match or fields_match.group(1).split("|") != CSQ_FIELDS:
        raise ValueError("Stage 2 CSQ header does not match the safety schema")
    counts = VCFCounts(header_file.header.contigs)

# Query VCF directly.
# No intermediate TSV is written.
fmt = (
    "%CHROM\t"
    "%POS\t"
    "%REF\t"
    "%ALT\t"
    "%INFO/CSQ\t"
    "%INFO/CLINVAR_ID\t"
    "%INFO/CLINVAR_CLNSIG\n"
)

proc = subprocess.Popen(
    [
        "bcftools",
        "query",
        "-f",
        fmt,
        INPUT_VCF,
    ],
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    text=True,
)


tmp_output = OUTPUT_PARQUET + ".tmp"

writer = pq.ParquetWriter(
    tmp_output,
    schema,
    compression="zstd",
)

buffer = []
batch_size = 50_000
input_variants = 0
output_rows = 0


def flush():
    global buffer, output_rows

    if not buffer:
        return

    table = pa.Table.from_pylist(
        buffer,
        schema=schema,
    )

    writer.write_table(table)

    output_rows += len(buffer)
    buffer = []


try:

    for line in proc.stdout:

        fields = line.rstrip("\n").split("\t")

        if len(fields) != 7:
            raise ValueError(
                f"Expected 7 columns, found {len(fields)}"
            )

        (
            chrom,
            pos,
            ref,
            alt,
            csq,
            clinvar_id_raw,
            clinvar_sig_raw,
        ) = fields

        counts.observe(chrom, pos, ref, alt)
        input_variants += 1

        records = parse_csq(csq)

        clinvar_id = clean(clinvar_id_raw)
        clinvar_sig = clean(clinvar_sig_raw)

        regulatory = regulatory_context(records)
        selected = selected_transcripts(records)

        # Intergenic / regulatory variant with no transcript
        if not selected:
            selected = [None]

        for record in selected:

            if record is None:
                gene_symbol = None
                gene_id = None
                transcript_id = None
                mane = None
                biotype = None
                consequence = None
                region = "Intergenic"

            else:
                gene_symbol = clean(record.get("SYMBOL"))
                gene_id = clean(record.get("Gene"))
                transcript_id = clean(record.get("Feature"))

                mane = (
                    "Yes"
                    if clean(record.get("MANE_SELECT"))
                    else None
                )

                biotype = clean(record.get("BIOTYPE"))
                consequence = clean(record.get("Consequence"))
                region = genomic_region(record)

            buffer.append({
                "Chr_Pos": f"{chrom}:{pos}",
                "Variant_REF_ALT": f"{ref}>{alt}",
                "Gene_symbol": gene_symbol,
                "Gene_ID": gene_id,
                "Transcript_ID": transcript_id,
                "MANE_Select": mane,
                "Transcript_biotype": biotype,
                "Consequence": consequence,
                "Genomic_region": region,
                "ClinVar_ID": clinvar_id,
                "ClinVar_significance": clinvar_sig,
                "Regulatory_context": regulatory,
            })

            if len(buffer) >= batch_size:
                flush()

    flush()

    stderr = proc.stderr.read()
    return_code = proc.wait()

    if return_code != 0:
        raise RuntimeError(
            "bcftools query failed:\n"
            + stderr
        )

finally:
    writer.close()


source_count = int(subprocess.check_output(["bcftools", "index", "-n", str(snakemake.input.source)], text=True).strip())
if source_count != input_variants or output_rows < input_variants:
    raise ValueError("Stage 2 annotation/expansion lost or duplicated variants")
record_qc("stage2_safety", input_variants, output_rows, input_variants, input_variants, "transcript expansion")

os.replace(
    tmp_output,
    OUTPUT_PARQUET
)


print(f"Input variants: {input_variants:,}")
print(f"Stage-2 transcript rows: {output_rows:,}")
print(f"Output: {OUTPUT_PARQUET}")
