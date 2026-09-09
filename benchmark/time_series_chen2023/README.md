# Barcode time series: Chen, Johnson, Hérissant et al. (2023)

The selected starting experiment is **hBFA1 in YPD**, a bulk fitness assay with
four sampled generations (8, 16, 24, 40) and two biological assay replicates.
The eight runs contain **16 paired FASTQ files, 1,248,297,814 compressed bytes
(1.25 GB)**. This subset is small enough for a controlled first cluster run.
It is a fitness-assay time series of pooled evolved clones, not the entire
long-term evolution experiment.

The [paper](https://doi.org/10.7554/eLife.92899) deposits raw sequencing under
[PRJNA912754](https://www.ebi.ac.uk/ena/browser/view/PRJNA912754).
The [author's analysis repository](https://github.com/mjohnson11/PLT_code/tree/a375a116bf69160f634a3d6a1b0ecc53ad62d142)
provides the sample/primer metadata and processed barcode counts. Its hBFA1
YPD reference contains 2,314 barcode pairs measured at the same eight samples.
None of these eight samples is excluded in the repository's timepoint-exclusion
table. The study also contains whole-genome and other assays; **do not download
the whole BioProject** for this first analysis.

The earlier Levy accession SRR5747458, used by Shepherd and Bartender, is a
single initial-timepoint run in PRJNA391509. It remains a useful clustering
benchmark but does not supply the requested deposited time series by itself.

## Selected samples

| Assay replicate | Generation | SRA run | Paired files, MB compressed |
|---|---:|---|---:|
| R1 | 8 | SRR22757105 | 105.9 |
| R1 | 16 | SRR22757108 | 148.8 |
| R1 | 24 | SRR22757107 | 264.7 |
| R1 | 40 | SRR22757106 | 138.4 |
| R2 | 8 | SRR22756607 | 135.6 |
| R2 | 16 | SRR22757104 | 158.3 |
| R2 | 24 | SRR22756609 | 145.0 |
| R2 | 40 | SRR22756608 | 151.6 |

Here R1/R2 in the **replicate** column are biological assay labels. FASTQ mate
1 and mate 2 are separate concepts. [samples.tsv](samples.tsv) records this
distinction, HTTPS URLs, ENA MD5 checksums and sizes, and extraction offsets.

## Read design and extraction

```text
FASTQ mate 1: UMI (8) — inline index (6) — fixed sequence — GGTACC — BC2 — ATAACT
FASTQ mate 2: UMI (8) — inline index (9) — fixed sequence — GGTACC — BC1 — ATAACT

BC2 = diverse lineage barcode; BC1 = environment/landing-pad barcode
Lineage identifier = BC2 + "_" + BC1; preserve the association between mates.
```

The expected component length is 26 bases, with variable blocks and internal
spacers. Both components are kept in their original read orientation, matching
the author's output. The downloaded pilot reads are **125 bases per mate**;
do not hardcode the paper's general 150-base sequencing description.

The selected FASTQs are already demultiplexed. The author's zero-based nominal
barcode starts are 63 in mate 1 and 49 in mate 2. The helper follows the
[published parser](https://github.com/mjohnson11/PLT_code/blob/a375a116bf69160f634a3d6a1b0ecc53ad62d142/processing/PLT_parse_predemultiplexed.py):

1. Check FASTQ integrity, equal mate counts, and matching read IDs.
2. Require mean Phred+33 quality ≥30 over the two nominal 26-base windows.
3. Search around each barcode for the author's flank patterns, permitting
   component lengths 24–28. Extract the observed sequence, retaining indels.
4. Deduplicate the concatenated 8+8-base UMI within each sample, retaining its
   first passing observation. Record duplicate UMIs and conflicting barcode
   pairs; this follows the paper's rule rather than using UMI consensus.
5. Write paired counts and separate component tables for barbac. Preserve all
   extracted pairs in the audit file; exclude pairs containing non-ACGT bases
   from the method-comparison inputs, recording how many molecules this removes.

The requested analysis uses **FastQC → PEAR → reference mapping → BAM → barbac
extraction**, with the existing `barbac_env`. The direct flank extractor below
is retained as a reproduction of the publication's parser and a validation
baseline. Its SLURM extraction stage currently runs that baseline; it does not
yet implement the requested BAM extraction of paired barcodes and UMIs.

The Chen cassette needs its own mapping reference. The older
`Reference_sequence_Barcdes` file in this project describes a different construct.
Mapping and extraction must preserve both barcode components, their association,
indels, and the sample's UMI rule before we compare clustering methods.

The [masked candidate FASTA](reference/chen2023_masked_amplicon.fasta) reconstructs
the 167-base cassette from modal fixed sequences in the pilot reads; it is not
an author-deposited reference. [Provenance](reference/provenance.json) records
the supporting counts and hashes. Nominal reference coordinates are BC2 50–75
and BC1 110–135 (one-based, inclusive); BC1 must be reverse-complemented for
the author's output orientation. These are nominal regions, not instructions
to discard insertions or force every extracted barcode to 26 bases.

The [complete-sample mapping pilot](mapping_pilot.md) ran the existing barbac
pipeline in the existing environment in about 4.1 minutes locally. It supports
local processing of the selected subset; it is not a full analysis runtime or
an extraction/clustering accuracy result.

The 24–28 interval deliberately matches the publication's extraction profile.
It does not retain arbitrarily large indels. A later broader extraction study
must use the same broader inputs for every clustering method and report its
extra yield separately.

## Local pilot completed

[Pilot results](pilot_result.json) use the first **50,000 complete read pairs**
of SRR22757105, obtained from small HTTP range downloads:

| Extraction outcome | Read pairs / molecules |
|---|---:|
| Input pairs | 50,000 |
| Failed quality threshold | 421 |
| Failed barcode pattern | 614 |
| Repeated UMI | 11 |
| Retained molecules | 48,954 |
| Distinct observed barcode pairs | 4,766 |

**97.908% is extraction retention, not clustering accuracy or F1.** The
retained counts agree exactly with the author's parser on the same pilot.
All 50,000 pairs carry the expected inline indices. Diverse barcode lengths
24, 25, 26, 27 and 28 are retained. The extraction took about 1.35 seconds
locally; this excludes downloading, FastQC and clustering, and is not a
cluster-runtime prediction. The prefix is not a random sample. No methods
have yet been compared on the complete time series.

## Existing environment and publication-parser reproduction on the cluster

Use the [project's existing environment instructions](../../README.md#environment-setup-optional--only-for-the-fastqbam-pipeline-native-install).
The barbac R package is installed in R; `configure_environment()` provisions
the external tools in `barbac_env`. A second `barbac-extract` environment is
unnecessary. The commands below currently reproduce the publication parser;
the final BAM extraction profile is still being validated.

The scripts are configurable because cluster partition/account and scratch
details have not yet been supplied. They use one CPU per sample, initially
8 GB RAM and two hours for FastQC plus extraction. Those are starting resource
requests, not measured requirements. Reserve about 10 GB of working space for
this subset and check `sacct` before expanding. No cluster jobs have been submitted.

1. Get the branch on the cluster and enter it:

   ```bash
   git clone --single-branch --branch feat/exact-search-clustering https://github.com/loukesio/barbac.git
   cd barbac
   ```

2. Install the working branch of barbac into the R library used on the cluster,
   if it is not already installed. Reuse the standard setup in R:

   ```r
   # Only if this branch is not already installed:
   remotes::install_github("loukesio/barbac@feat/exact-search-clustering")
   library(barbac)
   configure_environment()  # One-time creation; reuses an existing barbac_env
   use_barbac_env()         # Adds its tools to PATH in this R session
   check_barbac_tools()
   ```

   For a shell or SLURM job, activate the same external-tool environment with
   `conda activate barbac_env`, or load equivalent cluster modules. The batch
   job must also be able to find `Rscript` and the R library containing barbac
   when running the mapping/BAM workflow. Activating Conda alone does not install
   the barbac R package. The publication-parser helper uses Python 3's standard
   library; check `python3 --version` in the batch environment.

3. Copy [cluster_config.example.sh](cluster_config.example.sh) to your own
   cluster configuration file. Set the absolute repository and scratch paths,
   account and partition. `BARBAC_ENV_SETUP` may point to a shell file that
   activates the environment inside the batch job. Leave scheduler fields empty
   only when your cluster supplies working defaults. Record the actual environment:

   ```bash
   source /absolute/path/cluster_config.sh
   mkdir -p "$BARBAC_WORK"
   conda list --explicit > "$BARBAC_WORK/environment-explicit.txt"
   ```

4. Submit the first sample as a pilot. It downloads and MD5-verifies the full
   first pair of FASTQs (~106 MB), runs FastQC, then extracts the first 100,000
   pairs. This cluster pilot is distinct from the local 50,000-pair prefix check:

   ```bash
   bash benchmark/time_series_chen2023/submit.sh /absolute/path/cluster_config.sh pilot
   ```

5. Inspect the job state, logs and `extraction.json`. After the pilot completes,
   submit all eight samples for complete extraction:

   ```bash
   bash benchmark/time_series_chen2023/submit.sh /absolute/path/cluster_config.sh full
   ```

   The download array is followed by an extraction array with an `afterok`
   dependency. Pilot and full extraction outputs use different directories.
   [SLURM documents these array and dependency options](https://slurm.schedmd.com/sbatch.html).
   If compute nodes have no internet, use a permitted download partition or
   run the download stage on the cluster's transfer node first:

   ```bash
   for i in {0..7}; do
     python3 benchmark/time_series_chen2023/run_sample.py download \
       --index "$i" --work-dir "$BARBAC_WORK"
   done
   ```

   The download job skips files only after verifying size and MD5, so submitting
   the normal workflow after transfer does not download them again. Interrupted
   transfers remain `.part` files and are restarted on retry. A completed file
   with a wrong checksum is reported for inspection. Completed extractions are
   protected against accidental overwrite; use a new work directory for reruns.

6. Once jobs finish, collect resource usage and combine quality reports:

   ```bash
   sacct -j YOUR_EXTRACTION_JOB_ID --format=JobID,State,Elapsed,MaxRSS,ExitCode
   multiqc "$BARBAC_WORK/qc" -o "$BARBAC_WORK/multiqc"
   ```

## Outputs and the next comparison

Each sample has `extracted/full/SAMPLE/` containing:

- `barcode_pairs.csv`: both observed components and molecule counts, including
  any non-ACGT pairs for audit.
- `pairs_barbac.csv`: shared ACGT-only paired input for comparing methods.
- `diverse_barbac.csv` and `environment_barbac.csv`: `barcode,counts,barcode_length`
  tables accepted by `super_cluster2()`.
- `extraction.json`: count reconciliation, lengths, UMI conflicts, inline-index
  checks, extraction time, file hashes, settings and provenance.

For the later comparison, cluster the **two components separately**, preserve
their membership maps, and reconstruct BC2/BC1 pair counts for each time point.
Do not treat concatenated 52-base pairs as ordinary Hamming inputs; the native
Hamming implementation supports at most 32 bases. Do not collapse the lineage
identity to BC2 alone. See [comparison_plan.md](comparison_plan.md).

Published reference files and provenance can be retrieved without downloading
FASTQs:

```bash
python3 benchmark/time_series_chen2023/fetch_sources.py
python3 benchmark/time_series_chen2023/build_manifest.py
```

The reference subset is written to `generated/sources/published_YPD.csv`.
Downloads of author files are pinned to a repository revision and SHA-256.
ENA is a live service; if its metadata snapshot changes, the checksum check
stops and the changed metadata must be reviewed. The committed `samples.tsv`
already contains the selected download manifest and does not require refreshing.
Large source tables and sequencing files are ignored by Git.

To reproduce the local prefix check after fetching the sources, install
NumPy/pandas in a local analysis environment and run
`python3 benchmark/time_series_chen2023/reproduce_pilot.py`. This downloads
only the two recorded 2 MB prefixes, verifies their hashes, and checks every
extracted pair count against the pinned, reviewed author parser. It deliberately
does not reuse the full archive's MD5 to describe these partial downloads.

Validation: six extraction regression tests pass, including indel retention,
UMI conflicts, paired-file mismatch and truncated FASTQs. The local pilot also
matches the author's parser exactly. Bash syntax and submission argument checks
are local checks; the cluster environment and full FastQC/SLURM execution remain
to be validated on the user's cluster.
