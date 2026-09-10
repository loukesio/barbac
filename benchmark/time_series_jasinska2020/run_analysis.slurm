#!/bin/bash
#SBATCH --job-name=barbac-ecoli
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=barbac-ecoli-%j.log

set -euo pipefail
: "${BARBAC_PROJECT:?Set BARBAC_PROJECT to the checked-out barbac repository}"
cd -- "$BARBAC_PROJECT"

# Load your site's R/Quarto modules and existing barbac environment before
# submission, or add site-specific module commands here. No account or
# partition is guessed. Stage the checksum-verified FASTQs before this job.
srun Rscript benchmark/time_series_jasinska2020/build_release.R
srun Rscript benchmark/time_series_jasinska2020/process_samples.R
srun Rscript benchmark/time_series_jasinska2020/validate_processing.R
srun Rscript benchmark/time_series_jasinska2020/cluster_populations.R
srun Rscript benchmark/time_series_jasinska2020/benchmark_release.R "$SLURM_CPUS_PER_TASK"
srun Rscript benchmark/time_series_jasinska2020/validate_results.R
srun Rscript benchmark/time_series_jasinska2020/r_report/build_report.R
