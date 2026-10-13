# Rare Disease WGS Variant Prioritization — Analysis Methodology

This workflow prioritizes candidate variants from **proband-only whole-genome sequencing (WGS)** for rare-disease analysis. Starting from the normalized variant set, it combines genomic annotation, population and clinical evidence, functional predictions, gene–disease knowledge, phenotype matching, inheritance and literature evidence into a structured set of candidate hypotheses for review.

Annotation becomes progressively more detailed as the analysis advances, while the final candidate table integrates the relevant evidence into a single interpretable view.

## Workflow overview

```text
Normalized WGS variants
        ↓  Minimal genomic annotation + broad population routing

Broad retained variant set
        ↓  Population, clinical, transcript and regulatory annotation

Rare / clinically relevant variants
        ↓  Variant-level evidence + Gene/disease-level evidence

Fully annotated candidate hypotheses
        ↓  Inheritance organization + weighted scoring

Ranked candidate table for manual review and final clinical interpretation
```

## Outputs and scoring

The main output is an **information-rich ranked candidate table**, where each row represents a `variant × disease × inheritance` hypothesis.

Two supporting evidence layers are retained:

- **Variant evidence detail** — detailed evidence associated with individual variants.
- **Gene–disease detail** — gene-, disease-, phenotype- and literature-level evidence.

Each candidate receives a transparent **0–100 prioritization score** based on:

| Component | Weight |
|---|---:|
| Variant pathogenicity | **30%** |
| Mode of inheritance | **20%** |
| Gene / disease evidence | **25%** |
| Phenotype / case relevance | **15%** |
| Rarity | **10%** |

The component scores remain visible alongside the `Final_score`, allowing the ranking to be interpreted from the underlying evidence. The score is used to prioritize candidates for manual review and does not replace final clinical interpretation.

## Core principles

- Integrate complementary evidence rather than relying on a single predictor or database.
- Increase annotation depth as the candidate set becomes more focused.
- Treat missing evidence as unknown rather than benign.
- Reuse gene- and disease-level evidence across relevant candidates.
- Keep automated prioritization separate from final human interpretation.

Small summaries under `results/qc/` record each case’s stage input/output rows and unique variants. Checks follow each stage’s expected filtering, expansion or enrichment behavior; reference-case totals are not universal targets. Missing required fields, malformed records and unintended data loss stop execution.

Stage 1 requires unchanged annotation counts on canonical chromosomes and every contig represented in the VEP output. Offline VEP may omit entire unplaced (`Un_*`), unlocalized (`*_random`), or decoy (`*_decoy`) contigs with no annotation features. This exception is accepted only when VEP's unsupported-contig warnings account for every omitted record. The Stage 1 QC records those exclusions separately from common-frequency filtering; absent annotation does not establish biological irrelevance. Annotated noncanonical records remain subject to the normal routing rules.

Stage 2 remote gnomAD annotation allows four retries per chromosome after the initial attempt, with delays of 30, 60, 120 and 240 seconds. Retries repeat the remote header/query operation after a read failure; schema, count and scientific validation errors remain fatal. Successful chromosome QC records the attempts and retries used. Exhausting five attempts stops the workflow.

Remote header checks have a two-minute limit. Each remote chromosome query has a two-hour limit, configurable through `stage2.gnomad_remote_timeout_seconds` (default `7200`). On timeout the complete query process group is terminated and reaped before the existing retry policy applies. These limits apply to newly started workers; an existing worker retains the code it loaded at launch.

---

# Environment and Resource Preparation

Before configuring a case, create the two project-local Conda environments: `.conda/envs/rare-wgs` for resource preparation and the main workflow, and `.conda/envs/spliceai-env` for SpliceAI and FASTA indexing. Resource setup installs VEP 116.1 inside the main environment; Stage 4 uses the explicit project-local predictor prefix.

Run the one-action resource setup described in [sources_setup.md](sources_setup.md). The default catalog automatically imports checksum-verified GenCC and ClinGen exports from `config/resource_snapshots/` and downloads, processes and indexes the other local resources. `config/resource_catalog.json` defines the intended resource profile; `data/resources/resource_manifest.json` records actual readiness, provenance and identities.

`READY` means the required local resources passed their checks. Remote population, conservation, annotation and literature services are accessed during analysis; local readiness does not guarantee their availability or freeze their responses. Versions and scope limitations are documented in [resource_provenance.md](resource_provenance.md); preparation details are in [resource_processing_reference.md](resource_processing_reference.md).

For installation commands and the estimated 8-12-hour complete analysis runtime with eight cores, see [README.md](../README.md). Environment creation and initial resource setup take additional time. These instructions describe the supported Conda path.

---

# Case Context Setup

Before starting the analysis, the user provides basic information about the proband and the clinical case:

- **Sex** — used for sex-chromosome inheritance.
- **ORPHA / MONDO disease IDs** — when a known or suspected disease is available.
- **HPO terms** — patient phenotype terms, stored for phenotype-based prioritization in later stages.
- **Genes of interest** — manually provided suspected genes when available.
- **disease_names** — full name of the suspected disease/s if the diagnosis is clear. 
- **Broad disease domains** — optionlal layer that defines the disease broad clinic area`neurodevelopmental`, `neurological`, `skeletal`, `cardiac`, `renal`, `metabolic`, `hematologic`, `immune`, `ophthalmic`, `hearing`, `endocrine`, `gastrointestinal`, `dermatologic`, `craniofacial`, `multisystem`. These domains are used to guide literature searches when a specific disease name is unavailable. 
- **hom_max_override** — adjustable gnomAD homozygote-count limit for recessive candidates and hemizygote-count limit for X-linked candidates (default `2`). AD singleton candidates still require zero homozygotes; other filtering criteria are unchanged.
- **expected inheritance** - defines the inheritance model(s) considered plausible for the case. It accepts one or more predefined values such as `AD`, `AD_de_novo`, `AR_homozygous`, `AR_compound_heterozygous`, `X_linked`, `Mosaic`, and `Mitochondrial`. 

Known ORPHA diseases are linked to their associated genes. MONDO diseases are mapped through their corresponding ORPHA entries. These genes, together with manually provided genes, form the **case-gene list** used during Stage 3.

HPO terms are preserved but do **not** define case genes at this stage.

---

# Stage 1 — Minimal Genomic Annotation and Routing

## Goal
Create the minimum annotation required to remove clearly common variants before more detailed analysis.

## Fields
1. **Gene_symbol** — overlapping gene when available. Example: `BRCA1`
2. **Chr_Pos** — GRCh38 chromosome and normalized genomic position. Example: `chr17:43071077`
3. **Variant_REF_ALT** — normalized reference and alternate allele. Example: `C>T`, `A>AT`
4. **Consequence** — basic VEP consequence. Example: `missense_variant`, `intron_variant`
5. **gnomAD_cache_AF** — genome allele frequency available directly from the VEP cache. Example: `3.2e-6`

## Routing
Variants continue when:

```text
gnomAD_cache_AF ≤ 0.25
or
gnomAD_cache_AF is missing
```

The threshold is used only as a broad routing filter. Current population evidence is obtained in Stage 2.

## Main sources
- **Normalized VCF** — normalized coordinates and alleles
- **Ensembl VEP cache** — gene, consequence and cached gnomAD genome AF

---

# Stage 2 — Genomic, Population and Clinical Safety Layer

## Goal

Add current population, clinical, transcript, regulatory and patient-genotype information and preserve the result as a **transcript-level safety layer before inheritance filtering**.

The safety layer allows later rescue if the focused analysis gives no clear answer or only a partial explanation, including searches for:

* a second recessive allele;
* variants removed by population filtering;
* noncoding/regulatory variants without a clear gene assignment.

**Unit:** `variant × transcript`
**Output:** `results/tables/stage2_safety_annotated.parquet`

## Fields

1. **Chr_Pos** — normalized GRCh38 position. Example: `chr17:43071077`
2. **Variant_REF_ALT** — normalized REF and ALT. Example: `C>T`
3. **FILTER** — variant quality status from the input VCF. Example: `PASS`
4. **Gene_symbol** — associated gene. Example: `BRCA1`
5. **Gene_ID** — Ensembl gene ID. Example: `ENSG00000012048`
6. **Transcript_ID** — Ensembl transcript ID. Example: `ENST00000357654.9`
7. **MANE_Select** — MANE Select status. Example: `Yes`
8. **Transcript_biotype** — transcript type. Example: `protein_coding`
9. **Consequence** — transcript consequence. Example: `missense_variant`
10. **Genomic_region** — transcript-relative location. Example: `Exon 14`
11. **gnomAD_AF** — overall gnomAD AF. Example: `3.2e-6`
12. **gnomAD_grpmax_AF** — highest ancestry-group AF. Example: `8.1e-6`
13. **gnomAD_grpmax_population** — corresponding ancestry group. Example: `nfe`
14. **gnomAD_AC** — alternate-allele count. Example: `3`
15. **gnomAD_AN** — total called alleles. Example: `152340`
16. **gnomAD_het** — heterozygous carriers. Example: `3`
17. **gnomAD_hom** — homozygous alternate carriers. Example: `0`
18. **gnomAD_hemi** — hemizygous carriers where applicable. Example: `0`
19. **ClinVar_ID** — ClinVar Variation ID from the source VCF. Example: `12345`
20. **ClinVar_significance** — compact ClinVar classification. Example: `Pathogenic/Likely_pathogenic`
21. **Regulatory_context** — overlapping Ensembl/ENCODE regulatory elements. Example: `Ensembl: Enhancer | ENCODE: dELS`
22. **GT** — proband genotype. Example: `0/1`
23. **GQ** — confidence score for the proband genotype call. Example: `99`
24. **DP** — total read depth. Example: `42`
25. **AD_ref** — REF-supporting reads. Example: `23`
26. **AD_alt** — ALT-supporting reads. Example: `19`
27. **VAF** — alternate-allele fraction. Example: `0.45`

```text
VAF = AD_alt / (AD_ref + AD_alt)
```

## Main sources

* **Normalized VCF** — coordinates and genotype
* **Ensembl VEP / MANE** — transcript annotation
* **gnomAD v4.1 Joint** — population frequencies and counts
* **ClinVar** — clinical significance
* **Ensembl Regulatory / ENCODE** — regulatory context

HGVS and rsID annotations are added to the smaller candidate set: HGVS in Stage 3 and rsIDs in Stage 4. They are not generated in the Stage 2 safety table.

---

# Stage 3 — Inheritance Filtering and Variant Collapse

## Goal

Collapse Stage 2 transcript-level annotations into **one row per normalized variant**, then retain variants that remain plausible based on quality, functional relevance, population rarity, patient genotype, inheritance pattern and case information.

The resulting candidate set is used for deeper annotation and prioritization in the following stages.

---

## Filtering scheme / workflow

```text
One row per variant
        ↓
QC filtering
        ↓
Functional relevance
        ↓
Genotype + population rarity filtering
with a case-aware approach
        ↓
Inheritance organization
        │
        ├── AD
        ├── AR homozygous
        ├── AR compound heterozygous
        ├── X-linked
        ├── Mosaic
        ├── Mitochondrial
        └── Secondary
        ↓
Stage 3 candidates
```

### 1. QC filtering

Variants must:

- have **FILTER = PASS**
- have **GQ ≥20**
- have **DP ≥10**

Clinical or case-related evidence does not rescue variants that fail QC.

### 2. Functional filtering

Variants are retained when their annotation indicates a potentially relevant functional effect, including:

- coding and protein-altering consequences such as missense, truncating, frameshift and in-frame variants;
- splice-related variants or intronic variants close to a splice boundary;
- supported promoter or enhancer variants, including selected regulatory variants close to a relevant transcription start site;
- selected functional noncoding consequences such as mature miRNA variants.

ClinVar **Pathogenic / Likely Pathogenic (P/LP)** variants are retained after passing QC even when they do not meet the specific functional-group criteria above.

P/LP handling requires a main classification containing only Pathogenic, Likely pathogenic, or both. Mixed benign, uncertain, other or conflicting classifications do not qualify; the original evidence is preserved.

### 3. Genotype and population rarity filtering — case-aware approach

Variants are evaluated according to the proband genotype and the inheritance patterns they could support.

Population rarity is based on the highest available overall or ancestry-group gnomAD allele frequency. Missing population information is treated as **unknown**, rather than as evidence that the variant is common.

A **case gene** is a gene already supported by the clinical context provided before the analysis:

- a gene associated with the entered **ORPHA disease**;
- a gene associated with a **MONDO disease mapped through ORPHA**;
- or a **manually entered gene of interest**.

Variants in case genes receive broader retention thresholds, because there is already independent clinical reason to examine that gene more carefully.

This is especially important for possible **second hits**. If a P/LP variant is already present in a gene, another heterozygous variant in the same gene can be retained under the broader compound-heterozygous thresholds, allowing a possible recessive second allele to remain in the candidate set.

| Candidate type | Standard gene | Case gene |
|---|---|---|
| **AD singleton** | AF `≤0.0001`, `<30 het`, `0 hom` | AF `≤0.001`, `<50 het`, `0 hom` |
| **X-linked heterozygous** | AF `≤0.0001`, `<30 het`, `0 hom` | AF `≤0.001`, `<50 het`, `0 hom` |
| **AR homozygous** | AF `≤0.001`, `<500 het`, `≤2 hom` | Same |
| **AR compound candidate** | each allele: AF `≤0.001`, `<500 het`, `≤2 hom` | each allele: AF `≤0.01`, `≤2 hom`* |
| **X-linked hemizygous** | AF `≤0.001`, `<500 het`, `≤2 hemi` | Same |
| **Mosaic** | AF `≤0.001`, DP `≥20`, ALT reads `≥5`, VAF `10-35%` | Same |
| **Mitochondrial** | mtDNA AF `≤0.005`, heteroplasmy `≥10%` | Same |
| **ClinVar P/LP** | No rarity cutoff; must pass QC | Same |

\*The broader compound threshold also applies when the gene already contains a retained P/LP variant.

Heterozygous variants first enter a broader candidate pool. After distance filtering and compound-pair resolution, variants without a compound partner must pass the stricter singleton AD thresholds; P/LP retention remains independent.

### Transcript context and distance filtering

After genotype and population rarity filtering, transcript-specific context and distance annotations are added to all remaining variants.

VEP provides transcript-level HGVS and genomic context, and GENCODE is used to calculate splice and transcription-start-site distances for each relevant transcript.

Transcript-dependent information is retained in compact variant-level fields while preserving the transcript association.

Examples:

```text
Gene_context:
ENST1=exon 4/12 | ENST2=intron 3/8

Splice_context:
ENST1=+120 | ENST2=-7

TSS_context:
ENST1=-850 | ENST2=+4200
```

Distance is used for filtering only when proximity is part of the functional interpretation:

- intronic splice candidates require distance **<=50 bp** from a relevant splice boundary;
- selected regulatory variants such as dELS require distance **<=10 kb** from a relevant TSS.

Stage 4 nearby-splice routing also uses **<=50 bp**; the final workbook displays nearest splice distance only at **<50 bp**.

Distance does not remove variants retained because they are in a **case gene**, are **ClinVar P/LP**, or represent a possible **second hit**. Variants that are not being evaluated through a splice or regulatory route are also not filtered based on distance.

Regulatory relevance is retained explicitly so that later stages can use the Stage 3 regulatory decision directly.

### Compound-heterozygous candidates

A possible compound-heterozygous candidate contains at least **two heterozygous variants at different genomic loci in the same gene**. Different alternate alleles at one chromosome/position do not form independent compound hits.

The compound thresholds are intentionally broader than the final singleton AD thresholds, because a variant that is too common to be convincing as an isolated dominant candidate may still be relevant as part of a recessive pair.

Because the analysis is **proband-only**, these remain **possible compound-heterozygous pairs**; phase is not established.

A retained variant may belong to **more than one inheritance category** when several interpretations remain possible.

### Secondary findings

Variants in ACMG Secondary Findings genes are evaluated separately.

A reportable Secondary finding requires:

- a gene included in the ACMG Secondary Findings list;
- a ClinVar **Pathogenic / Likely Pathogenic (P/LP)** variant;
- a genotype compatible with the reported inheritance pattern.

Simply being located in an ACMG Secondary Findings gene is recorded as information and does not by itself retain or report the variant.

Secondary status may coexist with another inheritance category.

---

## Fields

1. **Case_gene_match** — indicates whether the variant is located in a gene associated with the known case or manually provided by the user.  
   Example: `True`

2. **Secondary_gene_match** — indicates whether the variant is located in an ACMG Secondary Findings gene. This alone does not mean that the finding is reportable.  
   Example: `True`

3. **Compound_het_partner** — the possible second variant when the candidate forms a potential compound-heterozygous pair.  
   Example: `chr2:166245732|G>A`

4. **inheritance_category** — one or more inheritance routes supported by the retained variant.  
   Examples:  
   `AD`  
   `AR_homozygous`  
   `AR_compound_heterozygous`  
   `X_linked_hemizygous`  
   `Mosaic`  
   `Mitochondrial`  
   `AD | Secondary`
5. **HGVSc** — transcript-specific coding HGVS.\
   Example: `ENST00000357654=c.5266dupC`

6. **HGVSp** — transcript-specific protein HGVS when applicable.\
   Example: `ENST00000357654=p.Gln1756ProfsTer74`

7. **Gene\_context** — transcript-specific genomic context.\
   Example: `ENST1=exon 14/23 | ENST2=intron 13/22`

8. **MANE_Transcript_ID** — existing MANE transcript identifier, retained from Stage 3 through scoring when available.\
   Example: `ENST00000349431.11`

9. **Splice\_context** — transcript-specific nearest splice-boundary distance.\
   Example: `ENST1=-7 | ENST2=+142`

10. **TSS\_context** — transcript-specific TSS distance.\
   Example: `ENST1=+850 | ENST2=-4200`

11. **Nearest\_splice\_distance** — shortest splice-boundary distance among the retained transcripts.\
    Example: `7`

12. **Nearest\_splice\_transcript** — transcript associated with `Nearest_splice_distance`.\
    Example: `ENST00000357654`

13. **Nearest\_TSS\_distance** — shortest TSS distance among the retained transcripts.\
    Example: `850`

14. **Nearest\_TSS\_transcript** — transcript associated with `Nearest_TSS_distance`.\
    Example: `ENST00000357654`

15. **Nearest\_TSS\_gene** — gene associated with the nearest retained transcript/TSS.\
    Example: `ENSG00000012048`

16. **Regulatory\_relevant** — indicates whether the variant satisfies the Stage 3 regulatory relevance logic, including any required distance criterion.\
    Example: `True`

17. **Compound_het_gene** — source gene or genes supporting the candidate compound relationships.
18. **Secondary_MOI** — inheritance information from the secondary-finding gene list.
19. **Secondary_reportable** — whether the candidate meets the supported secondary-finding classification and genotype rules.

Source compound genes and possible partners remain available even if Stage 6 selects a different gene hypothesis. A variant may belong to more than one category when several interpretations remain possible.

---

## Output

Stage 3 produces a **CSV candidate table with one row per retained variant**.

It contains the variants that remain plausible after QC, functional, population, genotype and inheritance filtering, together with the relevant information collected in the earlier stages.

This focused candidate set is then used for the more detailed variant-, gene-, phenotype- and disease-level analyses performed in the following stages.

The complete Stage 2 safety layer remains available separately for later rescue or targeted searches.

---

## Main sources

- **Stage 2 safety layer** — provides the variant, transcript, genotype, population, clinical and regulatory annotations entering Stage 3.
- **Ensembl VEP** — provides HGVSc, HGVSp and transcript-level genomic context; the same annotation run also preserves CADD, REVEL, SIFT, protein/domain, LOFTEE and NMD evidence for Stage 4.
- **GENCODE** — provides gene and transcript structure used for splice and transcription-start-site context.
- **gnomAD** — provides population allele frequencies and heterozygous, homozygous and hemizygous carrier counts.
- **ClinVar** — identifies known Pathogenic / Likely Pathogenic variants used for P/LP retention and Secondary Findings evaluation.
- **ACMG Secondary Findings resource** — defines the genes and inheritance patterns evaluated as potential Secondary Findings.
- **Orphanet** — links known ORPHA diseases, and MONDO diseases mapped through ORPHA, to their associated genes.
- **User-provided case context** — provides proband sex, known disease information and manually suspected genes used for case-aware filtering.

# Stage 4 — Variant-Level Evidence

## Goal

Add detailed **variant-level evidence** only for variants retained after Stage 3.

Some annotations are added to **all variants when available**, while others are calculated only for relevant variant types. A variant may belong to more than one annotation group, for example missense + splice.

The final output remains **one row per variant**, but transcript-dependent annotations retain the transcript association.

Missing evidence remains **unknown, not benign**.

## Variant-type routing

Variants are assigned to one or more evidence groups based on their transcript consequence and context:

- **Missense** — exonic missense variants.
- **LoF** — frameshift, stop-gained and essential splice loss-of-function variants.
- **In-frame indel** — exonic in-frame insertions or deletions.
- **Splice** — explicit splice-related consequences or variants within 50 bp of a relevant splice boundary.
- **Regulatory** — variants retained as regulatory-relevant in Stage 3.
- **Mitochondrial** — mitochondrial variants.

Different transcripts of the same variant may support different groups, and all relevant routes are retained.

Evidence sources are accessed only when the corresponding variant group is present.

## Variant-type annotation

| Variant type                             | Main evidence/tools                                    |
| ---------------------------------------- | ------------------------------------------------------ |
| All variants                             | dbSNP, ClinVar, CADD, phyloP100way, phastCons100way    |
| Missense                                 | REVEL, AlphaMissense, SIFT, protein/domain context     |
| Frameshift / stop-gained / LoF           | LOFTEE, NMD, protein context                           |
| In-frame insertion/deletion              | protein/domain context                                 |
| Splice / splice-region / nearby intronic | SpliceAI + splice junction context                                    |
| Coding variant near splice junction      | Coding evidence + SpliceAI                             |
| Regulatory / promoter / enhancer         | ReMM, regulatory context                               |
| Mitochondrial                            | Available common annotations; partial evaluation |

## Fields

1. **Known\_variant** — dbSNP rsID obtained in Stage 4 and displayed as `rsID` in the final workbook.\
   Example: `rs80357906`

2. **ClinVar\_conditions** — disease-specific ClinVar evidence.\
   Example: `Disease A = Pathogenic | Disease B = VUS`

3. **CADD** — general deleteriousness.\
   Example: `PHRED=28.6`

4. **phyloP100way** — evolutionary conservation.\
   Example: `5.4`

5. **phastCons100way** — conservation probability.\
   Example: `0.99`

6. **Protein\_context** — protein/domain context, retained per relevant transcript.\
   Example: `ENST1:Q1756P | BRCT domain`

7. **NMD\_context** — predicted NMD behavior, transcript-specific where applicable.\
   Example: `ENST1:Predicted NMD`

8. **Missense\_evidence** — missense predictors.\
   Example: `REVEL=0.87 | AlphaMissense=0.94 | SIFT=0.01`

9. **Splice\_prediction** — splice predictors.\
   Example: `SpliceAI=0.91`

10. **Regulatory\_prediction** — noncoding predictors.\
    Example: `ReMM=0.94`

11. **LoF\_evidence** — LoF/NMD evidence.\
    Example: `LOFTEE=HC | NMD=Predicted NMD`

12. **Mitochondrial\_evidence** — unavailable MT-specific evidence remains blank.

13. **is\_missense** — indicates that at least one retained transcript supports a missense route.\
    Example: `True`

14. **is\_lof** — indicates that at least one retained transcript supports a LoF route.\
    Example: `True`

15. **is\_inframe\_indel** — indicates an in-frame insertion/deletion route.\
    Example: `False`

16. **is\_splice\_relevant** — indicates an explicit or nearby splice-relevant route.\
    Example: `True`

17. **is\_regulatory** — indicates a regulatory-relevant route.\
    Example: `False`

18. **is\_mitochondrial** — indicates a mitochondrial variant.\
    Example: `False`

19. **Evidence\_routes** — all applicable evidence groups.\
    Example: `COMMON | MISSENSE | SPLICE`

20. **Transcript\_routes** — transcript-specific routes when they differ.\
    Example: `ENST1:MISSENSE | ENST2:SPLICE`

Mitochondrial candidates continue with available annotations when MT-specific evidence or transcript/gene context is missing. Unavailable evidence remains blank and is reported with a warning.

## Transcript-specific annotation rule

Transcript-dependent evidence is calculated independently for each relevant transcript and must not be collapsed by simply selecting a single value.

The final variant-level row may summarize multiple transcript annotations, but the association between **transcript → evidence** must remain recoverable.

MANE may be used to identify a preferred transcript for display, but it must not remove relevant evidence from other transcripts.

## Main sources

- **Ensembl VEP** — CADD, REVEL, SIFT, protein/domain context, LOFTEE and NMD.
- **dbSNP / NCBI Variation Services** — known variant identifiers such as `rsID`.
- **ClinVar** — condition-specific clinical significance and associated diseases.
- **UCSC phyloP / phastCons** — evolutionary conservation scores.
- **AlphaMissense** — missense pathogenicity score and prediction.
- **SpliceAI** — predicted splice gain/loss and splice-impact score.
- **ReMM** — pathogenicity/deleteriousness score for regulatory and noncoding variants.

MT-specific interpretation sources are not configured in the current workflow.

### Predictors not used in the current implementation

The current implementation intentionally does not use **Pangolin**, **Sei**, or **MutPred-Indel**.

- **Pangolin** is omitted because SpliceAI provides the primary splice-prediction evidence and the additional installation/runtime complexity was not justified for the current workflow.
- **Sei** is omitted because ReMM together with the retained regulatory and target-gene/TSS context provides the current regulatory evidence layer without requiring a separate heavy model environment.
- **MutPred-Indel** is omitted because its current integration and usage requirements are not suitable for a simple reproducible automated workflow. In-frame indels are instead evaluated using protein/domain context, sequence-change severity and gene/regional intolerance.

These predictors may be reconsidered in a future extended implementation. Their absence is intentional and must not be treated as missing variant-level evidence.

---

# Stage 5 — Gene, Disease, Phenotype and Literature Evidence

## Goal

Add a **gene-level evidence layer** that supports interpretation of gene–disease associations, gene constraint, phenotype relevance, inheritance/mechanism and literature evidence for the retained candidates.

When several genes occur in one row, every gene is analyzed independently and each field is reported as `GENE1: value || GENE2: value`.

Evidence-source errors are reported as warnings and leave the affected evidence blank; valid candidates continue through analysis.

## Fields

1. **Disease_associations** — diseases associated with the gene and the strongest available GenCC classification.  
   Example: `Bartter disease type 3 — Definitive`

2. **Missense_constraint** — gnomAD missense observed/expected ratio and Z-score.  
   Example: `O/E=0.893 | Z=1.36`

3. **LOEUF** — gene-level loss-of-function observed/expected upper-bound fraction; lower values indicate greater depletion of loss-of-function variants. 
   Example: `1.39`

4. **ClinVar_pathogenic_fraction** — P/LP fraction among recognised P/LP, benign/likely benign and VUS entries for the gene, counted by distinct ClinVar Variation ID. Gene-level evidence includes recognised low-penetrance P/LP labels and is separate from pure-main variant P/LP eligibility.  
   Example: `0.127628`

5. **ClinVar_pathogenic_count** — number of unique Pathogenic / Likely pathogenic ClinVar variants.  
   Example: `85`

6. **Known_variant_landscape** — Pathogenic / Likely pathogenic ClinVar variants grouped by molecular consequence.  
   Example: `P/LP: splice=15 | frameshift=26 | nonsense=18 | missense=23 | synonymous=1 | other=2`

7. **Gene_HPO_terms** — disease-specific HPO profiles associated with the gene.  
   Example: `Bartter syndrome type 3: HP:0001949 Hypokalemia; HP:0200114 Metabolic alkalosis`

8. **Matched_HPO_terms** — patient terms matched to gene-associated disease phenotypes by exact identity or an HPO ancestor/descendant relationship.
   Example: `HP:0001949 Hypokalemia <-> HP:0001949 Hypokalemia`

9. **Inheritance_and_mechanism** — inheritance patterns and mechanisms reported for the gene–disease associations.  
   Example: `Bartter disease type 3: Autosomal recessive`

10. **Literature_support** — publication counts by search layer, total reviewed publications and detected study types.  
   Example: `Gene-association publications retrieved: 5 | Total publications reviewed: 5 | Study types detected: human-cases=3, functional=1, animal=1`

11. **Top_publications** — up to three representative publications with PMID, PMCID or DOI, year and title.  
    Example: `PMID:41458276 | 2025 | Complexities of Bartter Syndrome Type III: A Case Study in Jordan.`

12. **Gene_symbol_transcript_aligned** — original transcript-aligned gene symbols preserved while the general `Gene_symbol` list is normalised and deduplicated. Downstream hypothesis selection and display use this alignment.

## Case context

Patient HPO terms are taken only from the terms manually provided in the Case Context.

The disease context used for literature searches is taken from `disease_names`. If no specific disease names are provided, `broad_disease_domains` are used instead. Specific disease names and broad domains are not combined.

Predefined case genes from the earlier ORPHA/MONDO mapping or manually supplied gene list remain marked as case genes.

If no patient HPO terms are available, no phenotype terms are inferred automatically from the disease or gene.

## Gene constraint

Gene sensitivity to variation is represented by:

- missense constraint from gnomAD;
- ClinVar pathogenic variant fraction and count;
- the known pathogenic variant landscape.

These measures are retained as separate evidence rather than combined into a single gene-level classification.

## Gene–disease evidence

A gene–disease association describes the evidence that variation in a gene is linked to a specific disease and mode of inheritance.

The strongest available GenCC classification is retained:

`Definitive > Strong > Moderate > Supportive > Limited`

Other possible classifications are `Disputed Evidence`, `Refuted Evidence`, `Animal Model Only` and `No Known Disease Relationship`.

Only `Moderate`, `Strong` and `Definitive` associations are used for literature searches.

## Literature evidence

Literature searches include:

1. **Exact variant** — gene + HGVSc, HGVSp or rsID.
2. **Variant–disease** — gene + variant + disease context.
3. **Gene–case context** — gene + specific disease or broad disease domain.
4. **Gene–disease association** — gene + a Moderate-or-stronger disease association.
5. **Rescue** — gene + mutation/variant + rare disease/genetic disorder.

All applicable specific searches are performed. Rescue is used only when they return no publications.

Disease-name alternatives may be used. Gene–disease terms must appear in the title or abstract, while exact-variant searches may also use article body text.

Results are deduplicated, cached and reused.

Responses must pass the hit-count, result-list and publication-record checks before being cached, including responses with zero matches. Failed or incomplete responses get three retries after the initial request. After exhaustion, successful search layers and usable partial records are retained and labeled `Partial evidence: search incomplete`; a failed review without usable records remains unavailable. Older unverified zero-result caches are rechecked when used. No retrieval dates are added to displayed evidence.

`Top_publications` prioritizes variant-specific, human-case, functional and animal evidence. Recent publications are preferred within each type; otherwise, Europe PMC relevance order is used.

## Main sources

- **GenCC** — curated gene–disease associations, source-specific classifications and inheritance information where available. PanelApp-derived submissions are accessed through GenCC rather than queried separately.

- **HPO** — disease phenotype profiles and matching to patient HPO terms.

- **ClinVar** — pathogenic variant fraction and count, known pathogenic variants and mutational landscape.

- **gnomAD constraint** — gene-level missense constraint.

- **Europe PMC** — literature search, structured literature summary and selection of key publications.

---

# Stage 6 — Weighted Candidate Scoring and Prioritization

## Goal

Rank candidates using a transparent **0–100 weighted score**.

Case-supported candidates additionally receive explicit review priority without modifying their underlying `Final_score`.

The final score is composed of five components:

| Component | Weight |
|---|---:|
| Variant pathogenicity | **30%** |
| Mode of inheritance | **20%** |
| Gene / disease evidence | **25%** |
| Phenotype / case relevance | **15%** |
| Rarity | **10%** |

### Scoring architecture 


Gene-level evidence is parsed once and reused across variants. For each variant:

1. Calculate `Variant_pathogenicity_score` for each relevant transcript.
2. Group transcripts by gene and retain the highest-scoring transcript for each gene.
3. Calculate the five Stage 6 component scores for each remaining gene candidate.
4. Calculate `Final_score` for each gene candidate.
5. Select the gene hypothesis using case-gene/pure-main P/LP priority, then Final Score, Variant Pathogenicity Score and Gene / Disease Evidence Score, with source gene order resolving otherwise equal hypotheses.
6. Write the winning scores back to the original variant row.

Evidence from different genes is never mixed. Stage 6 preserves exactly one output row per input variant.

#### Compound-heterozygous candidates

Pairs defined in Stage 3 are evaluated in addition to the single-variant candidates. For each gene, the strongest valid pair is selected, and its `Compound_pair_score` is the mean of the two candidate scores. If one member is assigned to a stronger alternative hypothesis, the next valid pair is considered; additional alleles remain compound candidates but do not receive a selected partner or pair score.

Case-gene and P/LP hypotheses retain priority during this comparison.

Same-gene additional alleles retain AR scoring even outside the selected main pair. AD fallback applies only when a different gene's AD hypothesis wins; inheritance, gene-disease compatibility, AD rarity and Final_score are recalculated, and selected compound fields are cleared. Absence of a reportable pair alone does not trigger AD fallback.

#### Scoring scales 

Continuous numerical predictors are converted to a `0–10` score using predefined numerical anchors and linear interpolation between them, then scaled to the component's assigned weight. Categorical evidence generally uses: 

- `10` = supporting 
- `5` = partial, uncertain or missing 
- `0` = not supporting / incompatible

#### Missing evidence

Missing applicable evidence generally contributes **50% of its possible points**, without counting as observed pathogenicity coverage. Important exceptions include no-gene gene/case scores of 0, case information without a supporting connection scoring 20, confirmed zero reviewed publications scoring 0, and unavailable MT rarity using the existing 50/100 default. See [Scoring_specification.md](Scoring_specification.md) for the complete rules.

## 1. Variant pathogenicity — 0–100

For nuclear variants:
- CADD: `20`
- Conservation: `10`
- Route-specific evidence: `70`

| Route | Variant-specific evidence |
|---|---|
| **Missense** | REVEL `30` + AlphaMissense `30` + protein context `10` |
| **Protein-altering in-frame indel** | functional-site evidence `50` + functional/domain context `20` |
| **Splicing** | SpliceAI `40` + splice context `30` |
| **LoF / truncating** | truncation/NMD context `40` + LOEUF `20` + MANE transcript relevance `10` |
| **Regulatory / noncoding** | regulatory annotation `30` + target-gene/TSS context `20` + ReMM `20` |

Each applicable route is scored independently and the highest route score is retained. 

Variant pathogenicity is calculated for each relevant transcript. Within each gene, the transcript with the highest complete pathogenicity score is retained. 

`Pathogenicity_evidence_coverage` records the fraction of applicable pathogenicity evidence that was observed.

**CADD** is converted to a 0–10 score and contributes up to 20 points: `≤5=0`, `10=2`, `15=4`, `20=6`, `22.5=8`, `≥25=10`. Values between anchors are linearly interpolated. Missing CADD receives `5/10`.

Mitochondrial candidates receive partial evaluation using available common pathogenicity evidence and the existing inheritance, gene/disease and case-relevance rules. A separate MT pathogenicity model is not implemented, and nuclear route-specific predictors are not substituted. The workbook identifies these candidates as partially evaluated.

### ClinVar adjustment

ClinVar applies an additive adjustment based on the main classification:

```text
Pure Pathogenic:                                    +20
Pure Likely pathogenic or pure P+LP:                 +15
P/LP mixed with non-P/LP or conflicting labels:        0
VUS without conflicting or other labels:              +5
Absent, other or conflicting alone:                    0
Benign:                                             -20
Likely benign or pure B+LB:                          -15
B/LB mixed with conflicting, other or VUS labels:    -10
```

Equivalent accepted spellings and separators use the same policy. Gene-level ClinVar summaries are evaluated separately; the complete adjustment rules are in [Scoring_specification.md](Scoring_specification.md).

The adjusted pathogenicity score is bounded to `0–100`.

The original ClinVar annotation remains available in the input table.

`Pathogenicity_evidence_coverage` records the fraction of applicable evidence observed.

## 2. Mode of inheritance — 0–100

Inheritance is scored against the inheritance model(s) expected for the case and defined in the Case Context.

- direct match → `100`
- matching compound-heterozygous configuration → `80`
- no expected inheritance information → `50`
- AD alternative selected from a compound candidate when AR is expected → `25`
- incompatible → `0`

A variant remains a compound candidate while it has a valid partner. It returns to an available single-variant route only when its compound hypothesis is not selected.

## 3. Gene / disease evidence — 0–100

Summarizes how strongly the candidate gene is established as a disease-associated gene and whether the observed candidate is compatible with the known disease mechanism.

Scoring components:

| Evidence | Points | 
|---|---:| 
| GenCC gene–disease strength | `50` | 
| MOI compatibility | `15` | 
| ClinVar gene landscape | `20` | 
| Literature support | `15` |

GenCC classifications include `Definitive`, `Strong`, `Moderate`, `Supportive` and `Limited`, as well as disputed, refuted, animal-only or no-known-relationship classifications. The strongest available GenCC gene–disease classification is selected and converted to numeric evidence.

MOI compatibility contributes up to 15 points: matching inheritance families receive `15`, incompatible families `0`, and missing/unrecognised information `7.5`. Mosaic candidates with AD gene evidence receive `7.5`. This component is separate from compatibility with the case's expected inheritance. 

ClinVar gene landscape combines: 
- matching P/LP molecular class (`10`); 
- total P/LP count (`5`); 
- pathogenic fraction (`5`). 

P/LP count scoring: `0=0`, `1=1`, `2–4=2`, `5–10=3.5`, `>10=5`, missing=`2.5`. 
Pathogenic fraction scoring: `0=0`, `>0–<0.05=1`, `0.05–<0.10=2.5`, `≥0.10=5`, missing=`2.5`. 

Literature support uses the number and type of relevant publications already collected in Stage 5. Human case reports provide the strongest support, followed by functional/animal evidence; absent literature scores `0`, while unavailable evidence is neutral.

## 4. Phenotype / case relevance — 0–100

`Case_relevance_score` measures how strongly a gene is connected to the specific case.

Evidence is evaluated hierarchically:

- predefined case gene → `100`
- direct disease context → `100`
- strong HPO overlap → high relevance
- relevant broad disease domain → intermediate relevance
- case information exists but no supporting connection is found → `20`
- no case-specific information is available → `50 / unknown`
- clearly unrelated context → low relevance

**Disease context** matches the case's disease names against the gene's disease associations, giving `100` for a match.

**Disease domain** matches the supplied broad domains against those associations, giving `70`. The higher applicable domain or HPO score is used. Literature contributes separately to gene/disease evidence.

When patient HPO terms are available, **HPO overlap** is the proportion of the gene's associated HPO terms that match the HPO terms provided for the case. General/ancestor HPO terms may also count as matches.

When no patient HPO terms are available, HPO evidence is treated as missing rather than negative.

Case relevance is gene-level evidence and can therefore be calculated once per gene and reused across variants.

## 5. Rarity — 0–100

Stage 3 defines the hard rarity thresholds; Stage 6 ranks the retained candidates according to how rare they are within those limits.

The rarity score combines:

- population AF support: `80%`
- inheritance-appropriate genotype-count support: `20%`

Population AF is scored relative to the applicable Stage 3 threshold:

- exactly at the threshold → `10`
- 5-fold below the threshold → `50`
- 50-fold below the threshold → `80`
- 500-fold or more below the threshold → `100`

Values between anchors are linearly interpolated. Zero or confirmed population absence receives `100`; missing AF receives `50`. When AF is missing but all relevant population genotype counts are explicitly zero, the variant is treated as absent from the population.

Genotype counts provide additional support according to the selected inheritance route. Zero counts provide the strongest support, high counts provide weaker support, and missing counts are neutral.

For compound-heterozygous candidates, each allele retains its own rarity score. Their complete candidate scores are combined in `Compound_pair_score`.

MT rarity is unavailable for scoring. The existing missing-evidence default gives `Rarity_score=50`, contributing five final-score points. MT frequencies remain available for review and are not treated as nuclear genotype counts.

## Final ranking

Candidates are grouped by `Scored_gene`, so variants assigned to the same gene remain together.

Groups containing a predefined case gene are displayed first. A gene with a selected compound pair is ranked using `Compound_pair_score`; otherwise, it is ranked using its highest `Final_score`.

Within each gene group, the selected compound pair is displayed first, followed by the remaining variants in descending `Final_score`.


## Fields added

1. **Variant_pathogenicity_score** — pathogenicity score of the transcript selected for the winning gene hypothesis. Example: `82`.
2. **Pathogenicity_evidence_coverage** — fraction of applicable pathogenicity evidence observed. Example: `0.75`.
3. **Inheritance_score** — compatibility with the expected inheritance model(s) of the case. 
4. **Gene_disease_evidence_score** — gene/disease-level evidence. Example: `80`
5. **Case_relevance_score** — relevance of the gene to the specific case. Example: `50`
6. **Rarity_score** — population rarity evidence. Example: `60`
7. **Scored_gene** — gene hypothesis that produced the retained `Final_score`. Example: `PBXIP1`
8. **Scored_transcript** — transcript hypothesis that produced the retained Variant_pathogenicity_score, including when Scored_gene is empty.Example: `ENST00000368463.8`
9. **Scored_route** — the variant route that was chosen for scoring in Variant_pathogenicity_score out of the possible variant routes. Example: `MISSENSE`
10. **`Scored_inheritance`** — inheritance model used for scoring.
11. **`AD_fallback`** — indicates that an AD hypothesis for a different gene won over the variant's AR compound-gene hypotheses. Additional alleles scored against the same AR gene retain AR inheritance, even when they are outside the selected main pair.
12. **`Scored_compound_partner`** — identifier of the other variant in the selected compound pair. Empty for AD fallbacks and variants without a selected pair.
13. **`Compound_pair_score`** — mean candidate score of the two variants in the selected compound pair. Empty for AD fallbacks and variants without a selected pair.
14. **`Final_score`** — weighted score of the selected per-variant hypothesis. Example: `88.10`.

---

# Stage 7 — Final Display and Human Review

## Goal

Convert the ranked analysis into a concise, human-readable review workbook without changing the biological evidence or scoring.

Stage 7 carries forward informative fields from Stages 2–6 and presents selected values in a compact format.

Transcript-level annotations and evidence are shown only for transcripts belonging to the `Scored_gene`. The `Scored_transcript` selected by Stage 6 is highlighted.

The selected transcript is blue; MANE transcripts are bold. `Gene_symbol` includes `| CASE` when the selected gene is listed in the case context.

## Review tabs and selection rules

The final workbook is `results/tables/rare_disease_variant_review.xlsx`. The default cutoff for review tabs is `Final_score ≥ 30`. The `P_LP` tab includes candidates whose main classification contains only Pathogenic, Likely pathogenic, or both, including those below the cutoff. Mixed or conflicting classifications are excluded from that tab. The `ALL` tab is unfiltered.

Tabs appear in this order:

* **Top_candidates** — includes case-gene candidates that meet the cutoff, then adds the highest-ranked non-case candidates matching the expected inheritance model. Each non-case variant must also have Variant_pathogenicity_score ≥ 30, including each allele of a compound pair. Case-gene variants are exempt from this additional pathogenicity minimum. The list aims to include at least 20 candidates, or 10 non-case candidates in addition to the qualifying case-gene candidates—whichever results in the larger list, subject to availability.
* **P_LP** — candidates meeting the main-classification P/LP rule above.
* **secondary** — reportable secondary findings.
* **AD** — autosomal dominant candidates.
* **AR_homozygous** — autosomal recessive homozygous candidates.
* **AR_compound** — valid compound-heterozygous pairs.
* **X_linked** — X-linked candidates.
* **Mosaic** — mosaic candidates.
* **Mitochondrial** — mitochondrial candidates.
* **ALL** — the complete candidate table.

The full workbook layout contains 42 columns; the included columns vary by tab.

## Compact field examples

Examples of fields that combine information for display:

* **Genotype** — `0/1 | 19/42 | 0.452` (GT, alternate reads/total depth, VAF).
* **Conservation** — `phyloP: 5.4 | phastCons: 0.99`.

Full field definitions are available in [field_reference_guide.csv](field_reference_guide.csv). Its column numbering follows the 42-column ALL layout.

## Human review

A human reviewer examines the candidates, assigns ratings as appropriate, and decides independently what is relevant. Automated scores guide review order; they do not determine the final biological conclusion.
