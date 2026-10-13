# Resource processing reference

These are maintenance details for `scripts/setup_resources.py`. For installation and recovery, use [sources_setup.md](sources_setup.md). The executable performs these steps automatically.

| Resource | Processing |
|---|---|
| Reference genome | Decompress, remove `chr` from sequence names, create FASTA `.fai` and pyfastx `.fxi` indexes. |
| ENCODE cCRE | Keep columns 1–3 and combine column 5 EH38E ID with column 6 class; sort chromosomes and positions; BGZF-compress and create a BED Tabix index. |
| ClinVar | Verify source MD5; keep chromosomes 1–22, X, Y and MT; rename MT to M; normalize and split alleles against the project FASTA; BGZF-compress and index. Reference allele mismatches cause failure. |
| AlphaMissense | Preserve canonical and isoform tables; create positional Tabix indexes. |
| gnomAD mitochondrial | Preserve the VCF and supplied index. |
| VEP cache | Extract safely to the configured cache directory; check assembly/GENCODE markers, chromosome folders and sample chunks; record cache identity; remove the archive after successful validation. |
| GENCODE, constraint, GenCC, HPO/HPOA, ClinGen, Orphadata and NMD plugin | Preserve source content and check the expected format. No biological filtering is applied. |
| ReMM index | Preserve the index and validate its configured checksum; score data remain remote. |
| SpliceAI models | Check the installed version and model files in the predictor environment. |

Downloads must pass configured checksums and basic format checks before use. Routine readiness checks inspect metadata, headers, required indexes and sample reads. They do not compare resources against another project or run analysis stages.

## Constraints when changing versions

- Keep the reference genome basename `GRCh38_no_alt_plus_hs38d1_no_chr.fasta`, the existing AlphaMissense filenames and `data/resources/orphanet/en_product1.xml` / `en_product6.xml`; some consumers discover these names directly.
- VEP software and cache major releases must match. Update `cache_version`, `assembly`, `cache_gencode` and the GENCODE `gtf_release_marker` consistently when changing those resources.
- ClinVar raw VCF, `.tbi` and `.md5` must come from the same snapshot. Keep Stage 5 raw-resource paths and the Stage 2 normalized path consistent with the catalog. Keep exactly one `clinvar_*.vcf.gz` in the active ClinVar directory; preserve earlier releases elsewhere.
- For ClinGen, update `report_marker` if the CSV title changes. All updated tables/XMLs must retain the fields expected by their consumers.
- For unchanged formats, a versioned `raw` filename may feed a stable `processed` path. Setup copies it automatically for resources of kind `copy`.
- Nuclear gnomAD uses the configured workflow URL template. Other remote URLs may also be embedded in consuming scripts; changing their catalog entry alone does not update those consumers.

Compatible configuration updates use the existing setup command. New schemas, assemblies, algorithms or services require a review of their consumers before use.
