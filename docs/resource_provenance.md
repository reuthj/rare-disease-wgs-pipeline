# Biological resource provenance

This is the resource version list for the default setup. Actual installed versions, readiness, download details and checksums are recorded in `data/resources/resource_manifest.json`. Source URLs and expected versions are defined in `config/resource_catalog.json`.

| Resource | Version or snapshot | Resource state |
|---|---|---|
| Reference genome | GRCh38 no-alt + hs38d1; GCA_000001405.15 analysis set | Local FASTA and indexes |
| GENCODE | v50 | Local GTF |
| ENCODE cCRE | Registry V4, GRCh38 | Local BED and index |
| VEP cache | Ensembl 116, GRCh38, GENCODE50 | Local extracted cache; archive removed after validation |
| gnomAD mitochondrial | Genomes v3.1 | Local VCF and index |
| ClinVar | GRCh38, 2026-09-13 | Local raw and normalized VCFs and indexes |
| AlphaMissense canonical | Zenodo record 8208688, hg38 | Local table and index |
| AlphaMissense isoforms | Zenodo record 8208688, hg38 | Local table and index |
| gnomAD constraint | v4.1.1 | Local table |
| GenCC | Preserved export; snapshot date unknown | Bundled pinned export; setup imports and verifies checksum |
| HPO ontology | 2026-09-01 | Local JSON |
| HPO disease annotations | 2026-09-02 annotations in tagged release v2026-09-01 | Local annotation table |
| ClinGen ACMG secondary findings | v3.3; export 2026-09-06 | Bundled pinned export; setup imports and verifies checksum |
| Orphadata product 6 | XML 2026-06-23; July 2026 knowledge base | Local disease–gene XML |
| Orphadata product 1 | XML 2026-06-23; July 2026 knowledge base | Local disease-mapping XML |
| NMD plugin | Ensembl release/116 | Optional local plugin |
| ReMM index | v0.4, hg38 | Local index for remote scores |
| SpliceAI models | 1.3.1 | Installed in the predictor environment |
| gnomAD nuclear | Joint v4.1 | Remote; accessed during analysis |
| ReMM scores | v0.4, hg38 | Remote; accessed during analysis |
| phyloP | hg38, 100way | Remote; accessed during analysis |
| phastCons | hg38, 100way | Remote; accessed during analysis |
| VEP REST | Live service; dataset release not pinned | Remote; accessed during analysis |
| NCBI variation / dbSNP | Live service; dbSNP release not pinned | Remote; accessed during analysis |
| Europe PMC | Live service; literature snapshot not pinned | Remote; accessed during analysis |

Local entries describe the prepared form, not proof that a particular installation is ready. The manifest records each installation's state. Remote service responses are not frozen by local resource setup.

MITOMASTER, MITOMAP, MitoTIP and HmtVar inputs for mitochondrial interpretation are not defined in this setup; their versions and readiness remain unspecified.
