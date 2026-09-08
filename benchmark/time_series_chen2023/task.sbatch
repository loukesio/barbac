#!/usr/bin/env bash
#SBATCH --job-name=barbac-extract
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
set -euo pipefail
config=$1
stage=$2
max_pairs=$3
source "$config"
if [[ -n ${BARBAC_ENV_SETUP:-} ]]; then source "$BARBAC_ENV_SETUP"; fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
kit="$BARBAC_REPO/benchmark/time_series_chen2023"
args=("$stage" --manifest "$kit/samples.tsv" --index "${SLURM_ARRAY_TASK_ID:?}"
      --work-dir "$BARBAC_WORK" --max-pairs "$max_pairs")
if [[ $stage == extract ]]; then args+=(--fastqc); fi
python3 "$kit/run_sample.py" "${args[@]}"
