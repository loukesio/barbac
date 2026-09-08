# Copy to a location on your cluster, edit, and pass that file to submit.sh.
# Paths must be absolute and visible on both login and compute nodes.
BARBAC_REPO=/path/to/barbac
BARBAC_WORK=/path/to/scratch/barbac_chen2023
BARBAC_ACCOUNT=""
BARBAC_PARTITION=""
BARBAC_DOWNLOAD_PARTITION=""
BARBAC_QOS=""
BARBAC_ARRAY_LIMIT=2
# Optional trusted shell file that activates Conda or loads Python/FastQC modules.
BARBAC_ENV_SETUP=""
# For example, that setup file could contain:
# source /path/to/miniforge3/etc/profile.d/conda.sh
# conda activate barbac-extract
