#!/bin/bash
# Submit the smoke test, then every experiment once it has passed.
# Usage (repo root, after 'make ARCH=-march=znver3'):  scripts/submit_all.sh [sbatch options, e.g. -p secondary]
set -e
mkdir -p results/logs
smoke=$(sbatch --parsable "$@" scripts/smoke.sbatch)
echo "smoke: $smoke"
for s in pingpong A_shape B_badP C_placement D_overlap E_packing; do
    echo "$s: $(sbatch --parsable -d afterok:"$smoke" "$@" scripts/$s.sbatch)"
done
