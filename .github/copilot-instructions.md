# Copilot instructions

## Project snapshot
- Investigates juvenile Pacific cod responses to thermal stress combining RNA-seq, WGS, methylation, and lipid assays.
- Workflows live in `code/` as numbered R Markdown notebooks that double as step-by-step lab logs.

## Repo layout & naming
- Preserve the `NN-description.Rmd` prefixes so downstream notebooks resolve expected outputs.
- Each notebook writes to `output/<same-name>/`; keep plots, logs, and rendered `.md`/`.html` in that folder.
- When adding notebooks, append a short entry to `code/README.md` pointing at its rendered output.

## Execution environment
- Heavy jobs assume UW Hyak paths (`/home/shared`, `/gscratch/scrubbed`); override variables in the `.bashvars` chunk before running elsewhere.
- Setup chunks generate `.bashvars` to centralize directories, program binaries, and thread counts—source it before executing bash steps.
- Most chunks set `eval = FALSE`; run commands manually in-shell to avoid firing long jobs during knitting.
- Keep program arrays (FastQC, MultiQC, Flexbar, kallisto, Trinity, ANGSD) intact so helper functions resolve the expected executables.

## Core pipelines
- `05-cod-RNAseq-trimming.Rmd` downloads reads from Owl, verifies MD5s, and stages FastQC/MultiQC plus Flexbar trimming.
- `06-cod-RNAseq-alignment.Rmd` handles kallisto indexing/quantification and links to HISAT2/featureCounts outputs; respect paths like `trimmed_reads_dir`, `kallisto_output_dir`.
- `07*` notebooks perform DESeq2 differential expression and GO analyses; they assume the kallisto matrices produced in `06`.
- `14-15` notebooks cover BWA alignments feeding ANGSD; BAM lists populate `output/15-angsd/bamlist.txt` automatically.
- `17.nextflow` runs nf-core/methylseq with `17.config` (Slurm queues, Singularity cache) and `17.samplesheet.csv`; outputs land in `output/17-nextflow`.

## Data management
- Raw reads and references live off-repo on Owl (`https://owl.fish.washington.edu/...`); download with the provided `wget --recursive --cut-dirs` pattern.
- Do not commit FASTQ/BAM; keep them on shared storage and reference via exported variables in `.bashvars`.
- Checksum every download with the MD5 snippets included in notebooks before proceeding downstream.

## Results & reporting
- Re-knit notebooks to refresh `.md`/`.html` after running commands, then push the updated prose alongside log summaries.
- Store MultiQC reports under the same output subtree (`raw-fastqc`, `trimmed-fastqc`) because later notebooks glob those directories.
- Document anomalies, reruns, and parameter tweaks directly in the R Markdown text for lab record continuity.

## Backups & housekeeping
- `bu.sh` rsyncs the repo (minus heavy artifacts and dotfiles) to Gannet; update exclusions if new large outputs appear.
- Singularity images cache at `/gscratch/scrubbed/srlab/.apptainer`; purge only after verifying no Nextflow rerun will need them.

## Validation expectations
- Treat successful MD5 checks, MultiQC summaries, and kallisto/DESeq2 logs as the “tests” to review before committing workflow changes.
- When tweaking upstream steps, regenerate dependent matrices and rerender the affected notebooks (e.g., `07-cod-RNAseq-DESeq2`) to confirm compatibility.
