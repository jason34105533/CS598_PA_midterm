#!/bin/bash
#SBATCH -J halo-E2
#SBATCH -p eng-instruction
#SBATCH -A 26fa-cs598pa-eng
#SBATCH --constraint=AE7713
#SBATCH -N 3
#SBATCH --ntasks-per-node=16
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH -t 00:30:00
#SBATCH -o results/logs/%x-%j.out
# Experiment E2 (bin/halo3d_ext): packing through index lists (indirect
# addressing), by hand (indexed) and as an MPI indexed-block type (indexed_dt).
# manual is rerun in the same job as the baseline.
cd "$SLURM_SUBMIT_DIR" && source scripts/common.sh
HALO_BIN=bin/halo3d_ext
[ -x $HALO_BIN ] || { echo "$HALO_BIN missing: run make"; exit 1; }

for N in 240 960; do
    for grid in 48x1x1 1x48x1 1x1x48 4x4x3; do
        for pack in manual indexed indexed_dt; do
            for rep in $(seq "$REPS"); do
                run 48 "$SPREAD" block --n "$N" --grid "$grid" --steps "$(steps_for "$N" 48)" \
                    --mode nonblocking --pack "$pack" --csv results/E2.csv --tag E2
            done
        done
    done
done
