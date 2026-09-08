#!/usr/bin/env bash
# Usage: bash submit.sh /absolute/path/cluster_config.sh pilot|full
set -euo pipefail
config=${1:?Supply an absolute path to your cluster config}
mode=${2:-pilot}
[[ $config == /* && -f $config ]] || { echo 'Config must be an existing absolute path' >&2; exit 1; }
source "$config"
[[ $BARBAC_REPO == /* && $BARBAC_WORK == /* ]] || { echo 'Use absolute repo/work paths' >&2; exit 1; }
kit="$BARBAC_REPO/benchmark/time_series_chen2023"
[[ -f $kit/samples.tsv ]] || { echo 'Sample manifest not found' >&2; exit 1; }
[[ ${BARBAC_ARRAY_LIMIT:-2} =~ ^[1-9][0-9]*$ ]] || { echo 'Invalid array limit' >&2; exit 1; }
case "$mode" in
  pilot) indices=0; max_pairs=100000 ;;
  full) indices="0-7%${BARBAC_ARRAY_LIMIT:-2}"; max_pairs=0 ;;
  *) echo 'Mode must be pilot or full' >&2; exit 1 ;;
esac
mkdir -p "$BARBAC_WORK/logs"
common=(--parsable --array="$indices" --chdir="$BARBAC_REPO"
        --output="$BARBAC_WORK/logs/%x-%A_%a.out" --error="$BARBAC_WORK/logs/%x-%A_%a.err")
[[ -z ${BARBAC_ACCOUNT:-} ]] || common+=(--account="$BARBAC_ACCOUNT")
[[ -z ${BARBAC_QOS:-} ]] || common+=(--qos="$BARBAC_QOS")
download_options=("${common[@]}")
extract_options=("${common[@]}")
partition=${BARBAC_DOWNLOAD_PARTITION:-${BARBAC_PARTITION:-}}
[[ -z $partition ]] || download_options+=(--partition="$partition")
[[ -z ${BARBAC_PARTITION:-} ]] || extract_options+=(--partition="$BARBAC_PARTITION")
download_job=$(sbatch "${download_options[@]}" --job-name=barbac-download --mem=1G --time=04:00:00 \
    "$kit/task.sbatch" "$config" download "$max_pairs")
download_job=${download_job%%;*}
extract_job=$(sbatch "${extract_options[@]}" --job-name=barbac-extract \
    --dependency="afterok:$download_job" "$kit/task.sbatch" "$config" extract "$max_pairs")
extract_job=${extract_job%%;*}
printf 'Download job: %s\nExtraction job: %s\nScope: %s\n' "$download_job" "$extract_job" "$mode"
printf 'Check usage after completion: sacct -j %s --format=JobID,State,Elapsed,MaxRSS,ExitCode\n' "$extract_job"
