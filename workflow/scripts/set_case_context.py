#!/usr/bin/env python3

from collections import defaultdict
from pathlib import Path
import xml.etree.ElementTree as ET

import yaml


# Isolated review output. Do not replace config/case_context.yaml until validated.
OUTPUT = Path("config/case_context.yaml")

ORPHA_GENE_FILE = Path("data/resources/orphanet/en_product6.xml")
MONDO_ORPHA_FILE = Path("data/resources/orphanet/en_product1.xml")

BROAD_DISEASE_DOMAINS = (
    "neurodevelopmental",
    "neurological",
    "skeletal",
    "cardiac",
    "renal",
    "metabolic",
    "hematologic",
    "immune",
    "ophthalmic",
    "hearing",
    "endocrine",
    "gastrointestinal",
    "dermatologic",
    "craniofacial",
    "multisystem",
)

EXPECTED_INHERITANCE_MODELS = (
    "AD",
    "AD_de_novo",
    "AR_homozygous",
    "AR_compound_heterozygous",
    "X_linked",
    "Mosaic",
    "Mitochondrial",
)

def ask_list(prompt):
    """Read a comma-separated list; empty input returns an empty list."""
    value = input(prompt).strip()
    if not value:
        return []
    return [x.strip() for x in value.split(",") if x.strip()]


def ask_broad_disease_domains():
    """Collect optional user-selected domains for focused literature searches."""
    print("\nAvailable broad disease domains:")
    print(", ".join(BROAD_DISEASE_DOMAINS))

    while True:
        values = ask_list(
            "\nBroad disease domains (comma-separated, optional):\n> "
        )
        if not values:
            return []

        normalized = [value.lower() for value in values]
        invalid = [
            value
            for value in normalized
            if value not in BROAD_DISEASE_DOMAINS
        ]
        if invalid:
            print("Invalid domain(s): " + ", ".join(invalid))
            print("Please use only values from the list above.")
            continue

        return list(dict.fromkeys(normalized))


def ask_sex():
    """Stage 3 currently requires known proband sex for X-chromosome routing."""
    while True:
        value = input("Proband sex [male/female]: ").strip().lower()
        aliases = {
            "m": "male",
            "male": "male",
            "f": "female",
            "female": "female",
        }
        if value in aliases:
            return aliases[value]
        print("Please enter: male or female.")

def ask_hom_max():
    """Ask for the recessive homozygote/X-hemizygote population count limit."""
    while True:
        value = input("Recessive hom/X-hemi count limit (Enter for default <=2; AD still requires 0 hom): ").strip()
        if not value:
            return None
        try:
            # Attempt to convert the input into an integer
            count = int(value)
            if count < 0:
                raise ValueError("Negative count")
            return count
        except (ValueError, TypeError):
            # Catch errors if it's text, a float, or an invalid type
            print("Enter a nonnegative integer, or press Enter for the default.")

def ask_expected_inheritance():
    """Collect one or more plausible inheritance models for the case."""
    print("\nAvailable inheritance models:")
    print(", ".join(EXPECTED_INHERITANCE_MODELS))

    while True:
        values = ask_list(
            "\nExpected inheritance models (comma-separated, optional):\n> "
        )

        if not values:
            return []

        aliases = {
            "ad": "AD",
            "ad_de_novo": "AD_de_novo",
            "de_novo": "AD_de_novo",
            "ar_homozygous": "AR_homozygous",
            "ar_hom": "AR_homozygous",
            "ar_compound_heterozygous": "AR_compound_heterozygous",
            "compound": "AR_compound_heterozygous",
            "compound_het": "AR_compound_heterozygous",
            "x_linked": "X_linked",
            "xl": "X_linked",
            "mosaic": "Mosaic",
            "mitochondrial": "Mitochondrial",
            "mito": "Mitochondrial",
        }

        normalized = []
        invalid = []

        for value in values:
            key = value.strip().lower().replace("-", "_").replace(" ", "_")
            model = aliases.get(key)

            if model is None:
                invalid.append(value)
            else:
                normalized.append(model)

        if invalid:
            print("Invalid inheritance model(s): " + ", ".join(invalid))
            print("Please use only values from the list above.")
            continue

        return list(dict.fromkeys(normalized))


def normalize_orpha_id(value):
    value = str(value).strip().upper().removeprefix("ORPHA:")
    if not value:
        return ""
    return f"ORPHA:{value}"


def normalize_mondo_id(value):
    value = str(value).strip().upper().removeprefix("MONDO:")
    if not value:
        return ""
    return f"MONDO:{value}"


def normalize_hpo_id(value):
    value = str(value).strip().upper().removeprefix("HP:")
    if not value:
        return ""
    return f"HP:{value}"


def find_orpha_code(disorder):
    """Handle capitalization variants seen across Orphadata products/versions."""
    for tag in ("OrphaCode", "ORPHAcode", "ORPHACode", "Orphacode"):
        code = disorder.findtext(tag)
        if code and code.strip():
            return code.strip()
    return None


def require_file(path):
    if not path.is_file():
        raise FileNotFoundError(f"Required resource not found: {path}")


def load_orpha_gene_map(path):
    """
    Build ORPHA ID -> associated gene symbols.

    Biological policy is intentionally unchanged from the original script:
    every gene association present in product6 is included. Association
    type/status are not filtered here.

    The parser supports the current GeneList layout and the older
    DisorderGeneAssociationList layout so a format change cannot silently
    turn the result into an empty map.
    """
    root = ET.parse(path).getroot()
    mapping = defaultdict(set)

    for disorder in root.findall(".//Disorder"):
        code = find_orpha_code(disorder)
        if not code:
            continue

        orpha_id = f"ORPHA:{code}"

        # Current Orphadata product6 layout: genes live under GeneList.
        gene_list = disorder.find("./GeneList")
        if gene_list is not None:
            for symbol_node in gene_list.findall(".//Symbol"):
                if symbol_node.text and symbol_node.text.strip():
                    mapping[orpha_id].add(symbol_node.text.strip())

        # Compatibility with older product6 layouts.
        for assoc in disorder.findall(
            "./DisorderGeneAssociationList/DisorderGeneAssociation"
        ):
            gene = assoc.findtext("./Gene/Symbol")
            if gene and gene.strip():
                mapping[orpha_id].add(gene.strip())

    if not mapping:
        raise ValueError(
            "No ORPHA-to-gene associations were parsed from en_product6.xml. "
            "The file may be wrong or its XML format may have changed."
        )

    return mapping


def load_mondo_orpha_map(path):
    """
    Build MONDO ID -> corresponding ORPHA IDs using Orphadata alignments.

    Biological policy is intentionally unchanged from the original script:
    all MONDO alignments present in product1 are accepted; mapping relation
    and validation status are not filtered here.
    """
    root = ET.parse(path).getroot()
    mapping = defaultdict(set)

    for disorder in root.findall(".//Disorder"):
        code = find_orpha_code(disorder)
        if not code:
            continue

        orpha_id = f"ORPHA:{code}"

        for ref in disorder.findall("./ExternalReferenceList/ExternalReference"):
            source = ref.findtext("Source")
            value = ref.findtext("Reference")

            if not source or not value:
                continue
            if source.strip().upper() != "MONDO":
                continue

            mondo_id = normalize_mondo_id(value)
            if mondo_id:
                mapping[mondo_id].add(orpha_id)

    if not mapping:
        raise ValueError(
            "No MONDO-to-ORPHA mappings were parsed from en_product1.xml. "
            "The file may be wrong or its XML format may have changed."
        )

    return mapping


def resolve_case_genes(
    orpha_ids,
    mondo_ids,
    manual_genes,
    orpha_gene_map,
    mondo_orpha_map,
):
    """
    Build Stage-3 case genes from known disease context and manual genes only.

    HPO terms are intentionally NOT converted to genes in Stage 3. They are
    retained in case_context.yaml for phenotype prioritization in later stages.
    """
    manual_set = {g.strip().upper() for g in manual_genes if g.strip()}
    case_genes = set(manual_set)

    source_genes = {
        "manual": set(manual_set),
        "orpha": set(),
        "mondo": set(),
    }
    unresolved = {
        "orpha": [],
        "mondo": [],
    }

    for orpha_id in orpha_ids:
        genes = set(orpha_gene_map.get(orpha_id, set()))
        if not genes:
            unresolved["orpha"].append(orpha_id)
        source_genes["orpha"].update(genes)
        case_genes.update(genes)

    for mondo_id in mondo_ids:
        mapped_orphas = set(mondo_orpha_map.get(mondo_id, set()))
        if not mapped_orphas:
            unresolved["mondo"].append(mondo_id)
            continue

        genes = set()
        for orpha_id in mapped_orphas:
            genes.update(orpha_gene_map.get(orpha_id, set()))

        if not genes:
            unresolved["mondo"].append(mondo_id)

        source_genes["mondo"].update(genes)
        case_genes.update(genes)

    return (
        sorted(case_genes),
        {source: len(genes) for source, genes in source_genes.items()},
        unresolved,
    )


def main():
    print("\n=== Stage 3 — Proband Context ===\n")
    print("Press Enter to skip any optional disease/phenotype field.\n")
    print("HPO terms are stored for later phenotype prioritization and do not create Stage-3 case genes.\n")

    for path in (ORPHA_GENE_FILE, MONDO_ORPHA_FILE):
        require_file(path)

    print("Loading local case-gene resources...")
    orpha_gene_map = load_orpha_gene_map(ORPHA_GENE_FILE)
    mondo_orpha_map = load_mondo_orpha_map(MONDO_ORPHA_FILE)

    print(f"ORPHA entries with genes:       {len(orpha_gene_map):,}")
    print(f"MONDO entries mapped to ORPHA:  {len(mondo_orpha_map):,}")

    sex = ask_sex()
    hom_max_override = ask_hom_max()
    expected_inheritance = ask_expected_inheritance()
    disease_names = ask_list(
        "\nDisease names (comma-separated, use full names if available, optional):\n> "
    )
    orpha_ids = sorted({
        normalize_orpha_id(x)
        for x in ask_list(
            "\nORPHA IDs (e.g. ORPHA:0001; comma-separated, optional):\n> "
        )
        if normalize_orpha_id(x)
    })
    mondo_ids = sorted({
        normalize_mondo_id(x)
        for x in ask_list(
            "\nMONDO IDs (e.g. MONDO:0000001; comma-separated, optional):\n> "
        )
        if normalize_mondo_id(x)
    })
    hpo_terms = sorted({
        normalize_hpo_id(x)
        for x in ask_list(
            "\nHPO terms (e.g. HP:0000001; comma-separated, optional):\n> "
        )
        if normalize_hpo_id(x)
    })
    broad_disease_domains = ask_broad_disease_domains()
    manual_genes = ask_list(
        "\nGenes of interest (gene symbols; comma-separated, optional):\n> "
    )

    case_genes, source_counts, unresolved = resolve_case_genes(
        orpha_ids,
        mondo_ids,
        manual_genes,
        orpha_gene_map,
        mondo_orpha_map,
    )

    print("\n=== Resolved case context ===")
    print(f"Sex: {sex}")
    print(f"disease names: {', '.join(disease_names) if disease_names else '-'}    ")
    print(f"ORPHA IDs: {', '.join(orpha_ids) if orpha_ids else '-'}")
    print(f"MONDO IDs: {', '.join(mondo_ids) if mondo_ids else '-'}")
    print(f"HPO terms: {', '.join(hpo_terms) if hpo_terms else '-'}")
    print(
        "Broad disease domains: "
        + (", ".join(broad_disease_domains) if broad_disease_domains else "-")
    )
    print(
        "Manual genes: "
        + (", ".join(sorted({g.upper() for g in manual_genes})) if manual_genes else "-")
    )

    print(
        "Expected inheritance: "
        + (", ".join(expected_inheritance) if expected_inheritance else "-")
    )

    print("\nGenes represented by each source (before de-duplication between sources)")
    print(f"Manual: {source_counts['manual']:,}")
    print(f"ORPHA:  {source_counts['orpha']:,}")
    print(f"MONDO:  {source_counts['mondo']:,}")
    print(f"Total unique case genes: {len(case_genes):,}")

    for source in ("orpha", "mondo"):
        if unresolved[source]:
            print(
                f"WARNING - unresolved {source.upper()} IDs: "
                + ", ".join(unresolved[source])
            )

    context = {
        "sex": sex,
        "hom_max_override": hom_max_override,
        "expected_inheritance": expected_inheritance,
        "disease_names": disease_names,
        "orpha_ids": orpha_ids,
        "mondo_ids": mondo_ids,
        "hpo_terms": hpo_terms,
        "broad_disease_domains": broad_disease_domains,
        "manual_genes": sorted({g.strip().upper() for g in manual_genes if g.strip()}),
        "case_genes": case_genes,
    }

    confirm = input(
        f"\nSave this configuration to {OUTPUT}? [Y/n]: "
    ).strip().lower()

    if confirm not in {"", "y", "yes"}:
        print("Cancelled.")
        return

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT, "w", encoding="utf-8") as f:
        yaml.safe_dump(
            context,
            f,
            sort_keys=False,
            allow_unicode=True,
        )

    print(f"\nSaved: {OUTPUT}")
    print(f"Final case genes: {len(case_genes):,}")


if __name__ == "__main__":
    main()
