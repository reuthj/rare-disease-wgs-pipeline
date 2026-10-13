# Biological resource setup

## Run setup

Use Conda on Linux x86_64 or WSL2. Create both environments once, then activate the main environment. Run these commands from the project root:

```bash
export CONDA_ENVS_PATH="$PWD/.conda/envs"
export CONDA_PKGS_DIRS="$PWD/.conda/pkgs"
conda env create --prefix "$PWD/.conda/envs/rare-wgs" -f envs/environment.yaml &&
conda env create --prefix "$PWD/.conda/envs/spliceai-env" -f envs/spliceai-env.yaml &&
conda activate "$PWD/.conda/envs/rare-wgs"
```

Run setup:

```bash
conda run --prefix "$PWD/.conda/envs/rare-wgs" --no-capture-output python scripts/setup_resources.py
```

Setup installs VEP 116.1 inside the main environment, then downloads, prepares and indexes local resources and reports **READY**. Run the same command again to resume an interrupted setup or repair incomplete resources. Valid resources are reused.

Plan for approximately **70 GB free** during setup and about **34 GB** for prepared resources. Initial resource preparation is estimated at **55–75 minutes**, depending on network and source availability. These are planning estimates. A complete Stage 1–7 whole-genome analysis with eight cores is estimated at **8–12 hours**, based on a roughly nine-hour run on the validation computer; environment creation and resource setup take additional time. Analysis outputs need additional space. As the final step, setup removes the VEP download archive after validating the extracted cache; its checksum and download details remain in the manifest.

Download connection attempts retry transient DNS/connection failures and HTTP 408, 429, 500, 502, 503 or 504 up to five times after the initial attempt (six attempts total), waiting 5, 10, 20, 40 and 60 seconds. Persistent failures still report INCOMPLETE; checksum and format checks remain mandatory.

Large downloads use **four parallel connections** by default. To change this:

```bash
conda run --prefix "$PWD/.conda/envs/rare-wgs" --no-capture-output python scripts/setup_resources.py --download-connections 8
```

Use `1` for sequential downloads; the supported range is 1–16. Servers without suitable support use one connection. This option is separate from Snakemake `--cores`.

The environment definitions provide the required tools, including the compiler used to install VEP. No separate VEP environment or system compiler installation is needed. `spliceai-env` provides SpliceAI and FASTA indexing. If your predictor environment has another name, use `--predictor-env NAME`. VEP software sources and revisions are recorded in `$CONDA_PREFIX/share/ensembl-vep-116.1/software_sources.json`. Existing compatible VEP installations are reused. Snakemake runs the same installer when creating rule environments with `--use-conda`.

## Check readiness

```bash
conda run --prefix "$PWD/.conda/envs/rare-wgs" --no-capture-output python scripts/setup_resources.py --check-only
```

This checks VEP and resources without downloading, installing, processing or removing archives.

| Result | Meaning |
|---|---|
| **READY** | All required local resources passed readiness checks. |
| **INCOMPLETE** | A required resource failed; setup exits with an error. |
| **SUBSET_READY** | Only resources requested with `--only` and their preparation dependencies passed. |

The generated `data/resources/resource_manifest.json` contains versions, status, checksums and failure details. Keep it with your analysis records. See [resource_provenance.md](resource_provenance.md) for the default version list.

Remote gnomAD, ReMM scores, conservation tracks and annotation/literature services are accessed during analysis. Their manifest state is `DEFERRED_TO_ANALYSIS_STAGE`; local READY does not check their availability. gnomAD data and cached remote indexes are organized under `data/resources/gnomad/` by resource and release. Case-specific chromosome annotations remain under `results/intermediate/stage2/gnomad_remote_by_chr/`. The local ReMM index is prepared by setup. Additional MITOMASTER/MITOMAP/MitoTIP/HmtVar inputs for mitochondrial interpretation are not covered by this setup.

## Recover an incomplete setup

Read the failed resource's `reason` and `next_step` in the manifest, correct the problem, then rerun the setup command.

| Problem | Action |
|---|---|
| Interrupted download or missing index | Rerun setup. Supported parallel transfers resume verified portions; sequential downloads restart. |
| Missing command or predictor environment | Create the environments above, activate `rare-wgs`, then rerun. |
| Interrupted VEP installation | Rerun setup; software is installed inside the active main environment. |
| Unreachable URL | Retry later or configure a mirror of the same file in the catalog. |
| Checksum mismatch | Supply the expected snapshot, or intentionally configure a compatible new version as described below. |
| Configuration path mismatch | Make the catalog and workflow configuration agree, then rerun. |

Do not delete the whole resource directory or remove checksums to bypass a failure. An unsuccessful VEP preparation keeps its archive for recovery. Old partial downloads or replaced incomplete cache directories may be retained with `.unverified-*` or `.incomplete-*` suffixes; remove them only after successful recovery.

The default GenCC and ClinGen exports are bundled in `config/resource_snapshots/`. Setup automatically validates the bundle and imports those exact exports; it downloads and prepares the other resources normally. Keep this directory with the project source. No separate resource preparation or snapshot argument is required for the default profile.

To intentionally supply another matching snapshot directory instead of the bundled files:

```bash
conda run --prefix "$PWD/.conda/envs/rare-wgs" --no-capture-output python scripts/setup_resources.py --snapshot-dir /path/to/snapshots
```

Files must match the catalog checksums. A damaged bundle fails setup rather than silently switching to a newer live export. Updated catalog profiles use their configured sources instead of the default bundle.

## Use a newer compatible version

Keep the original catalog, resource manifest and any snapshots needed to reproduce earlier results. For comparisons, use a separate project copy so the original resources remain available.

1. Copy `config/resource_catalog.json` to `config/updated_resources.json` and give its `profile` a new name.
2. Update the selected resource's `version`, `sources` and checksum together. Update release-specific fields when needed. A mirror of the same file needs only a URL change.
3. Keep `config/config.yaml` paths and configured URLs consistent with the new catalog. Keep the genome assembly, file format and required columns compatible with the workflow. VEP software, cache and GENCODE must use compatible releases. Changing VEP requires updating the pinned sources in `scripts/install_vep.py` as well as the configuration and resource catalog; changing only the cache URL is insufficient.
4. Run setup with the new catalog, then check it using the same catalog:

```bash
conda run --prefix "$PWD/.conda/envs/rare-wgs" --no-capture-output python scripts/setup_resources.py --catalog config/updated_resources.json
conda run --prefix "$PWD/.conda/envs/rare-wgs" --no-capture-output python scripts/setup_resources.py --catalog config/updated_resources.json --check-only
```

Use published checksums where available. If an intentionally updated resource has no configured checksum, setup records its first accepted file identity and rejects later content changes under the same specification. A different file format or remote service may require code changes; editing a URL alone does not establish compatibility.

See [resource_processing_reference.md](resource_processing_reference.md) for processing details and update constraints. New resource versions can change analysis results.
