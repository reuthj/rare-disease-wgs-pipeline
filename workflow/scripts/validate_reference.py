import re


header_file = snakemake.input.header
fai_file = snakemake.input.fai
output_file = snakemake.output[0]


# ---------------------------------------------------------
# Read contigs from VCF header
# ---------------------------------------------------------

vcf_contigs = {}

with open(header_file) as f:
    for line in f:
        if not line.startswith("##contig=<"):
            continue

        id_match = re.search(
            r"ID=([^,>]+)",
            line
        )

        length_match = re.search(
            r"length=([0-9]+)",
            line,
            flags=re.IGNORECASE
        )

        if not id_match:
            continue

        contig = id_match.group(1)

        length = (
            int(length_match.group(1))
            if length_match
            else None
        )

        vcf_contigs[contig] = length


# ---------------------------------------------------------
# Read contigs from FASTA index
# ---------------------------------------------------------

ref_contigs = {}

with open(fai_file) as f:
    for line in f:
        fields = line.rstrip().split("\t")

        contig = fields[0]
        length = int(fields[1])

        ref_contigs[contig] = length


# ---------------------------------------------------------
# Compare
# ---------------------------------------------------------

# Every contig used by the VCF must exist in the reference
# with exactly the same name.
missing_in_reference = (
    set(vcf_contigs)
    - set(ref_contigs)
)

# Extra contigs in the reference are allowed.
extra_in_reference = (
    set(ref_contigs)
    - set(vcf_contigs)
)

length_mismatches = []

for contig in (
    set(vcf_contigs)
    & set(ref_contigs)
):

    vcf_length = vcf_contigs[contig]
    ref_length = ref_contigs[contig]

    # Some VCF headers may omit contig length.
    # In that case existence is checked,
    # but length cannot be validated.
    if (
        vcf_length is not None
        and vcf_length != ref_length
    ):
        length_mismatches.append(contig)


compatible = (
    not missing_in_reference
    and not length_mismatches
)


# ---------------------------------------------------------
# Write report
# ---------------------------------------------------------

with open(output_file, "w") as f:

    f.write("metric\tvalue\n")

    f.write(
        f"vcf_contigs\t"
        f"{len(vcf_contigs)}\n"
    )

    f.write(
        f"reference_contigs\t"
        f"{len(ref_contigs)}\n"
    )

    f.write(
        f"missing_in_reference\t"
        f"{len(missing_in_reference)}\n"
    )

    f.write(
        f"extra_in_reference\t"
        f"{len(extra_in_reference)}\n"
    )

    f.write(
        f"length_mismatches\t"
        f"{len(length_mismatches)}\n"
    )

    f.write(
        f"compatible\t"
        f"{compatible}\n"
    )

    # Useful diagnostic details
    if missing_in_reference:
        f.write(
            "missing_contigs\t"
            + ",".join(
                sorted(missing_in_reference)
            )
            + "\n"
        )

    if length_mismatches:
        f.write(
            "length_mismatch_contigs\t"
            + ",".join(
                sorted(length_mismatches)
            )
            + "\n"
        )


if not compatible:
    raise ValueError(
        "Reference FASTA is not compatible with the VCF. "
        "Contig names and lengths must match. "
        "See reference validation report."
    )