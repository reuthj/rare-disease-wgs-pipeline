import os
import re


# ---------------------------------------------------------
# Inputs / outputs supplied automatically by Snakemake
# ---------------------------------------------------------

vcf_file = snakemake.input.vcf
header_file = snakemake.input.header
samples_file = snakemake.input.samples
stats_file = snakemake.input.stats

summary_file = snakemake.output.summary
fields_file = snakemake.output.fields


# ---------------------------------------------------------
# Helper functions
# ---------------------------------------------------------

def extract_id(line):
    """
    Extract ID from VCF header definitions such as:

    ##INFO=<ID=DP,...>
    ##FORMAT=<ID=GT,...>
    """

    match = re.search(r"ID=([^,>]+)", line)

    if match:
        return match.group(1)

    return None


def detect_genome_build(reference, contigs):
    """
    Try to identify the genome build.

    Header information is preferred.
    Chromosome lengths are used only as supporting evidence.
    """

    reference_lower = reference.lower()

    # Stronger evidence: explicit reference name
    if "grch38" in reference_lower or "hg38" in reference_lower:
        return "GRCh38", "reference header"

    if "grch37" in reference_lower or "hg19" in reference_lower:
        return "GRCh37", "reference header"

    # Weaker evidence: chromosome 1 length
    chr1_length = contigs.get("chr1", contigs.get("1"))

    if chr1_length == 248956422:
        return "GRCh38-like", "chromosome 1 length"

    if chr1_length == 249250621:
        return "GRCh37-like", "chromosome 1 length"

    return "Unknown", "no reliable evidence"


# ---------------------------------------------------------
# Read VCF header
# ---------------------------------------------------------

with open(header_file) as handle:
    header_lines = [line.rstrip() for line in handle]


fileformat = "Unknown"
reference = "Not specified"

sources = []

contigs = {}

info_ids = []
format_ids = []
filter_ids = []

header_fields = []


for line in header_lines:

    # VCF version
    if line.startswith("##fileformat="):
        fileformat = line.split("=", 1)[1]

    # Reference genome
    elif line.startswith("##reference="):
        reference = line.split("=", 1)[1]

    # Variant caller / software source
    elif line.startswith("##source="):
        sources.append(
            line.split("=", 1)[1]
        )

    # Contigs
    elif line.startswith("##contig=<"):

        contig_id_match = re.search(
            r"ID=([^,>]+)",
            line
        )

        length_match = re.search(
            r"length=([0-9]+)",
            line,
            flags=re.IGNORECASE
        )

        if contig_id_match:

            contig_id = contig_id_match.group(1)

            if length_match:
                contig_length = int(
                    length_match.group(1)
                )
            else:
                contig_length = None

            contigs[contig_id] = contig_length

    # INFO definitions
    elif line.startswith("##INFO=<"):

        field_id = extract_id(line)

        if field_id:

            info_ids.append(field_id)

            header_fields.append(
                ("INFO", field_id, line)
            )

    # FORMAT definitions
    elif line.startswith("##FORMAT=<"):

        field_id = extract_id(line)

        if field_id:

            format_ids.append(field_id)

            header_fields.append(
                ("FORMAT", field_id, line)
            )

    # FILTER definitions
    elif line.startswith("##FILTER=<"):

        field_id = extract_id(line)

        if field_id:

            filter_ids.append(field_id)

            header_fields.append(
                ("FILTER", field_id, line)
            )


# ---------------------------------------------------------
# Read sample information
# ---------------------------------------------------------

with open(samples_file) as handle:

    samples = [
        line.strip()
        for line in handle
        if line.strip()
    ]


# ---------------------------------------------------------
# Parse bcftools stats
# ---------------------------------------------------------

stats = {}


with open(stats_file) as handle:

    for line in handle:

        # SN = Summary Numbers in bcftools stats
        if not line.startswith("SN"):
            continue

        fields = line.rstrip().split("\t")

        if len(fields) >= 4:

            metric = fields[2].rstrip(":")
            value = fields[3]

            stats[metric] = value


# ---------------------------------------------------------
# Determine chromosome naming style
# ---------------------------------------------------------

if contigs:

    chr_prefixed = sum(
        contig.startswith("chr")
        for contig in contigs
    )

    if chr_prefixed > len(contigs) / 2:
        contig_style = "chr-prefixed"
    else:
        contig_style = "non-chr-prefixed"

else:
    contig_style = "Unknown"


# ---------------------------------------------------------
# Genome build inference
# ---------------------------------------------------------

genome_build, build_evidence = detect_genome_build(
    reference,
    contigs
)


# ---------------------------------------------------------
# Check whether annotation already exists
# ---------------------------------------------------------

has_vep = "CSQ" in info_ids
has_snpeff = "ANN" in info_ids

clinvar_fields = [
    field
    for field in info_ids
    if field.startswith("CLN")
]


# ---------------------------------------------------------
# Important genotype fields
# ---------------------------------------------------------

important_format_fields = [
    "GT",
    "DP",
    "AD",
    "GQ",
    "PL"
]

available_genotype_fields = [
    field
    for field in important_format_fields
    if field in format_ids
]


# ---------------------------------------------------------
# File metadata
# ---------------------------------------------------------

file_size_mib = (
    os.path.getsize(vcf_file)
    / (1024 ** 2)
)


# ---------------------------------------------------------
# Construct summary
# ---------------------------------------------------------

summary = {

    "input_vcf":
        vcf_file,

    "file_size_MiB":
        f"{file_size_mib:.2f}",

    "vcf_format":
        fileformat,

    "reference_header":
        reference,

    "inferred_genome_build":
        genome_build,

    "build_evidence":
        build_evidence,

    "contig_style":
        contig_style,

    "number_of_contigs":
        len(contigs),

    "number_of_samples":
        len(samples),

    "sample_names":
        ",".join(samples),

    "analysis_mode":
        snakemake.config["analysis"]["mode"],

    "source":
        ";".join(sources)
        if sources
        else "Not specified",

    "number_of_INFO_fields":
        len(info_ids),

    "number_of_FORMAT_fields":
        len(format_ids),

    "number_of_FILTER_fields":
        len(filter_ids),

    "important_genotype_fields":
        ",".join(available_genotype_fields),

    "has_VEP_CSQ":
        has_vep,

    "has_SnpEff_ANN":
        has_snpeff,

    "ClinVar_like_fields":
        ",".join(clinvar_fields)
        if clinvar_fields
        else "None",
}


# ---------------------------------------------------------
# Add selected bcftools statistics
# ---------------------------------------------------------

interesting_stats = [

    "number of records",
    "number of SNPs",
    "number of indels",
    "number of MNPs",
    "number of multiallelic sites",
]


for metric in interesting_stats:

    if metric in stats:
        summary[metric] = stats[metric]


# ---------------------------------------------------------
# Write main summary
# ---------------------------------------------------------

with open(summary_file, "w") as handle:

    handle.write("metric\tvalue\n")

    for metric, value in summary.items():

        handle.write(
            f"{metric}\t{value}\n"
        )


# ---------------------------------------------------------
# Write detailed header-field table
# ---------------------------------------------------------

with open(fields_file, "w") as handle:

    handle.write(
        "category\tfield_id\tdefinition\n"
    )

    for category, field_id, definition in header_fields:

        handle.write(
            f"{category}\t{field_id}\t{definition}\n"
        )
