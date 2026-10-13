# 🧬 Rare Disease WGS Variant Prioritization

![Genome](https://img.shields.io/badge/Genome-GRCh38-2ea44f)
![Workflow](https://img.shields.io/badge/Workflow-Snakemake-f39c12)
![Annotation](https://img.shields.io/badge/Annotation-Ensembl%20VEP-7b61ff)
![Language](https://img.shields.io/badge/Python-3.x-3776AB)
![Input](https://img.shields.io/badge/Input-Proband%20VCF-lightgrey)

A reproducible Snakemake workflow for prioritizing candidate variants from **proband-only whole-genome sequencing (WGS)** in rare-disease analysis,
starting from an already generated **VCF file**

The workflow integrates genomic and transcript annotation, population frequency, functional predictions, clinical evidence, gene–disease associations, inheritance models, phenotype information and literature evidence into a structured pipeline that progressively reduces a whole-genome variant set to a ranked collection of candidate hypotheses for manual review.

> **Core principle:** preserve evidence first, prioritize later.


## 🧭 Navigation

- [✨ Why this workflow?](#-why-this-workflow)
- [🧬 Project workflow](#-project-workflow)
- [📊 Scoring model](#scoring-model)
- [📥 Inputs & case configuration](#-inputs--case-configuration)
- [📚 Required resources](#-required-resources)
- [📤 Candidate outputs](#-candidate-outputs)
- [✅ Quality control & reproducibility](#-quality-control--reproducibility)
- [📁 Project structure](#-project-structure)
- [⚙️ Installation](#️-installation)
- [▶️ Running the workflow](#️-running-the-workflow)
- [⚠️ Limitations](#️-limitations)
- [📄 Documentation](#-documentation)


## ✨ Why this workflow?

The project is designed to function both as a rare-disease variant prioritization workflow and as a reusable, reproducible bioinformatics pipeline:

- **Reproducible** — software dependencies, workflow logic, thresholds and reference-resource versions are explicitly defined and tracked.
- **Modular** — analysis stages are separated into Snakemake steps, allowing individual parts of the workflow to be rerun without unnecessarily repeating completed upstream analyses.
- **Evidence-preserving** — available evidence is retained throughout the pipeline while filtering and prioritization progressively reduce the candidate space.
- **Hypothesis-aware** — variants are evaluated in the context of gene, transcript and inheritance hypotheses rather than by variant pathogenicity alone.
- **Transparent** — filtering decisions and prioritization scores are based on explicit rules and evidence components that remain available for inspection.
- **Scalable** — annotation, filtering, evidence integration and scoring are separated so individual resources or methods can be updated without redesigning the complete workflow.
- **Privacy-aware** — patient-specific data and generated analysis results remain local and are separated from the reusable workflow code.

---

## 🧬 Project workflow

### Pipeline overview

```text
🧬 Input VCF
   │
   ▼
🧹 Prepare & normalize
   │
   ▼
🧩 Stage 1 — Minimal Genomic Annotation and Routing
   │
   ▼
🛡️ Stage 2 — Base Annotation & Safety Layer
   │
   ▼
🧬 Stage 3 — Filtering & Inheritance Routing
   │
   ▼
🔬 Stage 4 — Variant-Level Evidence
   │
   ▼
📚 Stage 5 — Gene & Disease Evidence
   │
   ▼
📊 Stage 6 — Integrated Candidate Scoring
   │
   ▼
🔎 Stage 7 — Final Report & Human Review
```
---

### What each stage does

| Stage | Purpose | Main output |
|---|---|---|
| **Preparation** | Inspect the VCF, validate the reference and normalize variants | Normalized VCF + QC |
| **1 — Minimal Genomic Annotation and Routing** | Add the minimum annotation required to remove clearly common variants and route the remaining variants before more detailed analysis | Early-filtered and routed VCF |
| **2 — Base Annotation & Safety Layer** | Add core genomic, genotype and population information while maintaining broad variant retention | Evidence-enriched safety table |
| **3 — Filtering & Inheritance Routing** | Apply QC, functional and rarity criteria and assign plausible inheritance hypotheses | Retained candidate variants |
| **4 — Variant-Level Evidence** | Add transcript-aware functional and pathogenicity evidence appropriate to each variant route | Variant evidence table |
| **5 — Gene & Disease Evidence** | Add gene-level disease, phenotype, constraint, clinical and literature context | Gene- and disease-enriched candidates |
| **6 — Integrated Candidate Scoring** | Score candidate hypotheses using explicit and interpretable evidence components | Ranked candidate table |
| **7 — Final Report & Human Review** | Organize ranked candidates and supporting evidence for efficient manual interpretation | Final candidate review workbook |

---

<a id="scoring-model"></a>
## 📊 Scoring model

After evidence integration, each candidate hypothesis receives a transparent prioritization score based on five evidence components:

Variant evidence       30%
Gene–disease evidence  25%
Inheritance fit        20%
Case relevance         15%
Rarity                 10%

The score provides a common scale for comparing candidate hypotheses that may be supported by very different types of evidence.

| Component | Weight | What it represents |
|---|---:|---|
| **Variant evidence** | **30%** | Functional and clinical evidence supporting a potentially damaging effect of the variant, evaluated according to the relevant variant route (e.g., missense, loss-of-function, in-frame indel, splice or regulatory)  |
| **Gene–disease evidence** | **25%** | Strength of the evidence linking the scored gene to human disease, including established disease associations and supporting gene-level evidence |
| **Inheritance fit** | **20%** | Compatibility between the candidate inheritance hypothesis and the expected or documented mode of inheritance |
| **Case relevance** | **15%** | Case-specific evidence supporting the gene as relevant to the investigated phenotype or prior case knowledge |
| **Rarity** | **10%** | Population rarity of the variant in the context of the candidate inheritance model |

The final score is used to prioritize candidate hypotheses for manual review rather than to classify variants as pathogenic or benign.

Each component remains available separately in the output, allowing the evidence driving the final ranking to be inspected directly.

Detailed scoring rules, thresholds, route-specific predictors, missing-evidence handling and score calculations are documented separately in the **Scoring Specification**.

---

## 📥 Inputs & case configuration

The workflow requires two case-specific inputs: the proband VCF and the case context.

### Input VCF

The analysis starts from a **single-sample proband VCF** generated using the **GRCh38 genome build**.

### Case context

Case-specific information is provided through:

```text
config/case_context.yaml
```

The case context defines the clinical and analytical information used for case-specific prioritization, including phenotype information, known case genes and the expected inheritance model when available.

After resource setup, run the interactive case configuration:

```bash
python workflow/scripts/set_case_context.py
```

It collects the case information, resolves case genes and writes `config/case_context.yaml` for the analysis.

---

## 📚 Required resources

The workflow uses local annotation datasets and remote services. Downloaded resources are stored locally and are not included in the repository.

### Main resources

The main resources used by the workflow include:

- **GRCh38 reference genome**
- **Ensembl VEP cache and plugins**
- **GENCODE**
- **gnomAD** population frequency and gene constraint data
- **ClinVar**
- **ENCODE cCRE regulatory annotations**
- **CADD**
- **REVEL**
- **AlphaMissense**
- **SpliceAI**
- **ReMM**
- **phyloP / phastCons conservation scores**
- **GenCC**
- **PanelApp-derived gene–disease evidence**
- **Human Phenotype Ontology (HPO)** and phenotype annotations
- **Orphadata** disease–gene associations and disease mappings
- **ClinGen ACMG secondary-finding gene list**

Resource setup and recovery instructions are provided in `docs/sources_setup.md`. Resource versions and snapshot states are listed in `docs/resource_provenance.md`.

To use compatible newer releases, follow the setup guide to update versions, URLs and checksums in `config/resource_catalog.json`, keeping the workflow configuration consistent.

---

## 📤 Candidate outputs

The final output is `results/tables/rare_disease_variant_review.xlsx`, an Excel workbook that organizes prioritized candidates into review tabs while preserving their supporting evidence.

### Workbook tabs

| Tab | Content |
|---|---|
| **Top_candidates** | Main review shortlist. Includes case-gene candidates and the highest-scoring additional candidates compatible with the expected inheritance model |
| **P_LP** | Candidates with a main classification of pathogenic and/or likely pathogenic |
| **secondary** | Reportable secondary-finding candidates |
| **AD** | Autosomal-dominant candidates passing the final reporting threshold |
| **AR_homozygous** | Homozygous autosomal-recessive candidates passing the final reporting threshold |
| **AR_compound** | Retained compound-heterozygous candidates, organized by gene and scored variant pair |
| **X_linked** | Retained X-chromosome candidates, including relevant heterozygous and hemizygous variants |
| **Mosaic** | Retained mosaic candidates |
| **Mitochondrial** | Reserved for mitochondrial candidates |
| **ALL** | Complete scored candidate set, sorted by final score |

Most inheritance-specific tabs use **Final_score >=30**; `AR_compound` uses a selected reciprocal pair score of **>=30**. `Top_candidates`, `P_LP` and `secondary` use their own selection rules.

### Main output fields

Each candidate row preserves the evidence required to understand why the candidate was retained and how it was scored. The main field groups include:

- **Variant identity** — chromosome, genomic position, REF/ALT and rsID.
- **Gene and transcript context** — scored gene symbol, transcript, transcript biotype, consequence, HGVS, gene context and regulatory or protein context.
- **Genotype and population evidence** — genotype, read support, variant allele fraction, gnomAD allele frequency and population genotype counts.
- **Inheritance hypothesis** — the inheritance category under which the candidate is evaluated, such as AD, AR homozygous, AR compound heterozygous, X-linked or mosaic.
- **Variant-level evidence** — ClinVar, CADD, conservation, functional evidence route and route-specific scores. Routes include missense, loss-of-function, in-frame indel, splice and regulatory evidence.
- **Gene-level evidence** — LOEUF, missense constraint, disease associations, known inheritance and disease mechanisms, ClinVar gene-level landscape and phenotype information.
- **Literature evidence** — summarized literature support and selected publications relevant to the candidate gene or disease.
- **Scoring fields** — variant pathogenicity, inheritance fit, gene–disease evidence, case relevance and rarity component scores, followed by the weighted `Final_score`.
- **Compound-heterozygous fields** — retained pairs additionally include the scored partner and pair-level score used for ranking the compound hypothesis.

For example, a candidate may contain:

```text
Gene_symbol:                  CBS
Inheritance_category:         AR_homozygous
Evidence_routes:              REGULATORY
gnomAD_AF:                    0
Variant_pathogenicity_score:  20.00
Inheritance_score:            100.00
Gene_disease_evidence_score:  100.00
Case_relevance_score:         20.00
Rarity_score:                 100.00
Final_score:                  64.00
```

The complete definition and interpretation of every output column is provided in the **Field Reference Guide**.

The automated ranking determines review order and preserves the evidence supporting each candidate. Final biological or clinical interpretation remains manual.

---

## ✅ Quality control & reproducibility

The workflow includes quality-control and consistency checks throughout the analysis to ensure that filtering, annotation and scoring behave as expected.

Key principles include:

- **Input validation** — the input VCF and reference genome are checked before downstream analysis.
- **Stage-level tracking** — variant and row counts are recorded across major filtering and processing steps to identify unexpected losses or expansions.
- **Evidence preservation** — supporting annotations remain available in downstream outputs so filtering and scoring decisions can be traced back to their underlying evidence.
- **Transcript consistency** — transcript-dependent annotations are kept aligned with the transcript hypothesis being evaluated.
- **Explicit missing data** — unavailable evidence is handled according to defined rules rather than silently interpreted as negative evidence.
- **Deterministic scoring** — prioritization is based on explicit thresholds, weights and rules, allowing the same inputs and resource versions to reproduce the same ranking.
- **Progress preservation & caching** — completed outputs and reusable intermediate results are preserved, allowing interrupted runs to resume without repeating successfully completed work.
- **Versioned resources** — external databases, reference files and annotation resources are tracked by version wherever applicable.
- **Workflow reproducibility** — Snakemake manages dependencies between stages and reruns only the steps affected by changed inputs or outputs.

---


## 📁 Project structure

```text
rare-disease-wgs/
├── config/
│   ├── config.yaml
│   ├── resource_catalog.json
│   └── resource_snapshots/
│       ├── manifest.json
│       └── preserved_gene_evidence_snapshots.tar.gz
├── data/
│   ├── processed/
│   ├── raw/
│   ├── reference/
│   └── resources/
├── docs/
│   ├── field_reference_guide.csv
│   ├── project_specification.md
│   ├── Scoring_specification.md
│   ├── resource_provenance.md
│   ├── resource_processing_reference.md
│   └── sources_setup.md
├── logs/
├── results/
│   ├── intermediate/
│   ├── qc/
│   ├── tables/
│   │   └── rare_disease_variant_review.xlsx
│   └── variants/
├── workflow/
│   ├── Snakefile
│   ├── rules/
│   └── scripts/
├── envs/
│   ├── environment.yaml
│   ├── environment.post-deploy.sh
│   └── spliceai-env.yaml
├── scripts/
│   ├── install_vep.py
│   └── setup_resources.py
├── .gitignore
└── README.md
```

### What each folder contains

| Folder | Contents |
|---|---|
| `config/` | Workflow settings, resource catalog and case-context example |
| `config/resource_snapshots/` | Pinned GenCC and ClinGen exports imported automatically during setup |
| `data/` | Local case inputs and prepared resources |
| `data/raw/` | User-provided proband VCF and index |
| `data/processed/` | Normalized input VCF and index |
| `data/reference/` | Reference genome, annotation files and VEP cache |
| `data/resources/` | Biological datasets and the resource readiness manifest |
| `docs/` | Workflow, scoring, output-field and resource guides |
| `logs/` | Execution logs used to investigate failures |
| `results/` | Generated analysis outputs |
| `results/intermediate/` | Intermediate evidence, checkpoints and reusable caches |
| `results/qc/` | Input checks and stage-level QC summaries |
| `results/tables/` | Candidate tables and the final review workbook |
| `results/variants/` | Retained variant files and candidate keys |
| `workflow/` | Main Snakemake entry point and analysis stages |
| `workflow/rules/` | Snakemake rules and stage dependencies |
| `workflow/scripts/` | Case setup, annotation, filtering, scoring and report scripts |
| `envs/` | Main and SpliceAI environment definitions and VEP installation hook |
| `scripts/` | Automated resource preparation and VEP installation |

The tree includes runtime directories created locally during setup and analysis. Patient inputs, downloaded resources, logs and generated results remain outside version control. Case configuration creates the local `config/case_context.yaml` from user-provided context.

---

## ⚙️ Installation

Use Conda on Linux x86_64 or WSL2. Create both environments and activate the main environment:

```bash
export CONDA_ENVS_PATH="$PWD/.conda/envs"
export CONDA_PKGS_DIRS="$PWD/.conda/pkgs"
conda env create --prefix "$PWD/.conda/envs/rare-wgs" -f envs/environment.yaml &&
conda env create --prefix "$PWD/.conda/envs/spliceai-env" -f envs/spliceai-env.yaml &&
conda activate "$PWD/.conda/envs/rare-wgs"
```

Run these commands from the project root. The environments and package cache are created locally under `.conda/`.

SpliceAI uses the exact project-relative prefix `.conda/envs/spliceai-env` in `stage4.spliceai_conda_prefix`. Setup also selects the predictor beside the running project-local main environment by exact prefix.

Run resource setup from the project root after activating the environment:

```bash
conda run --prefix "$PWD/.conda/envs/rare-wgs" --no-capture-output python scripts/setup_resources.py
```

Setup downloads, processes and indexes resources, then reports **READY**. Remote services are accessed during analysis.

The project includes a small, checksum-verified bundle containing the exact GenCC and ClinGen exports used by the default catalog. Setup imports these automatically from `config/resource_snapshots/` and downloads/prepares the other resources. No separate snapshot download or `--snapshot-dir` argument is needed.

Plan for approximately **70 GB of free disk space** for resource setup, including download archives, and about **34 GB** for prepared resources. Initial resource preparation is estimated at **55–75 minutes**, depending on network and source availability. These are planning estimates; analysis outputs require additional disk space.

As a final step, setup removes the VEP download archive after validating the extracted cache. Its checksum and download provenance remain in the resource manifest.

Downloads use **four parallel connections** by default. Adjust this with `--download-connections N`; use `1` for sequential downloading.

To check readiness without downloading or processing:

```bash
conda run --prefix "$PWD/.conda/envs/rare-wgs" --no-capture-output python scripts/setup_resources.py --check-only
```

If setup reports **INCOMPLETE**, inspect the failed resource in `data/resources/resource_manifest.json`, follow the recovery instructions in `docs/sources_setup.md`, and rerun setup. Valid resources are reused. Archived exports can be supplied with `--snapshot-dir /path/to/snapshots`.

Review the provided workflow configuration and adjust only the local settings required for your environment, such as resource or input paths:

```text
config/config.yaml
```

Place your VCF (with index csi file) in `data/raw/` and set `input_vcf` in `config/config.yaml` to its actual filename. For example:

```yaml
input_vcf: "data/raw/your_sample.vcf.gz"
```

Case-specific information is configured separately using:

```bash
python workflow/scripts/set_case_context.py
```

Once the environment and required resources are prepared and the necessary local paths are set, the workflow is ready to run.

---

## ▶️ Running the workflow

The workflow is executed through Snakemake. Dependencies between stages are resolved automatically, so only the steps required to produce the requested output are executed.

### Dry run

Before starting an analysis, inspect the planned jobs without executing them:

```bash
conda run --prefix "$PWD/.conda/envs/rare-wgs" --no-capture-output snakemake --snakefile workflow/Snakefile --cores 8 -n
```

### Full run

Run the complete workflow:

Keep the computer plugged in and awake, with sleep and hibernation disabled and internet access available throughout the run.

```bash
conda run --prefix "$PWD/.conda/envs/rare-wgs" --no-capture-output snakemake --snakefile workflow/Snakefile --cores 8
```

For a complete Stage 1–7 whole-genome analysis with eight cores, allow approximately **8–12 hours**, based on a roughly nine-hour run on the validation computer. Runtime varies with hardware, input size, network and remote-service availability. Environment creation and initial resource setup take additional time.

The number of cores can be adjusted according to the available system resources.

Completed outputs are reused automatically, allowing the workflow to continue from existing results rather than repeating successfully completed stages.

### Run to a specific target

Individual workflow targets can also be requested directly:

```bash
conda run --prefix "$PWD/.conda/envs/rare-wgs" --no-capture-output snakemake --snakefile workflow/Snakefile --cores 8 <target>
```

Snakemake will execute only the requested target and any required upstream dependencies.

---

## ⚠️ Limitations

The current workflow has several important limitations:

- **Proband-only analysis** — the workflow is currently designed for single-proband analysis and does not perform trio- or family-based segregation analysis.
- **VCF-level starting point** — the workflow starts from an already generated VCF and does not include raw-read processing, alignment or variant calling.
- **Limited mitochondrial analysis** — dedicated mitochondrial variant calling and mitochondrial-specific analysis tools are not currently included.
- **No CNV or structural variant analysis** — the current workflow focuses on small variants and does not analyze copy-number variants (CNVs) or structural variants (SVs).
- **Large external resource dependency** — the workflow relies on multiple genomic, population, clinical and functional annotation resources, making the initial setup relatively complex.
- **Long runtime** — whole-genome annotation and evidence integration can require substantial processing time, particularly during the first complete run.
- **High storage requirements** — local reference genomes, annotation databases, caches and intermediate files require substantial disk space. Allow approximately **70 GB** of free space for initial resource setup, plus space for analysis outputs.

These limitations describe the current implementation and may be addressed in future versions of the workflow.

---

## 📚 Documentation

Additional documentation is available in the `docs/` directory:

- **`project_specification.md`** — detailed specification of the workflow stages, analysis logic and evidence collected throughout the pipeline.
- **`Scoring_specification.md`** — complete definition of the scoring system, including component weights, thresholds, route-specific scoring and handling of missing evidence.
- **`field_reference_guide.csv`** — reference guide for the fields included in the final candidate report, including their meaning, interpretation and relevant value ranges.
- **`sources_setup.md`** — resource setup, compatible version updates, archived snapshots and failure recovery.
- **`resource_provenance.md`** — resource versions, snapshots and local or remote state.

Together, these documents provide the detailed technical and analytical reference for reproducing, reviewing and extending the workflow.
