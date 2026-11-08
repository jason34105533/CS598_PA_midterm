#!/bin/bash
#SBATCH -J halo-C
#SBATCH -p eng-instruction
#SBATCH -N 4
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH -t 01:00:00
#SBATCH -o results/logs/%x-%j.out
# Experiment C: same grid, different rank placement.
#   block_spread : consecutive ranks on one node, one rank per CCD (baseline)
#   cyclic_spread: ranks dealt round-robin over nodes, one rank per CCD
#   block_packed : consecutive ranks on one node, 16 ranks packed on 2 CCDs
cd "$SLURM_SUBMIT_DIR" && source scripts/common.sh

for N in 256 1024; do
    for grid in 4x4x4 1x1x64 64x1x1; do
        for rep in $(seq "$REPS"); do
            common=(--n "$N" --grid "$grid" --steps "$(steps_for "$N")" --mode nonblocking --pack manual --csv results/C.csv)
            run 64 "$SPREAD" block "${common[@]}" --tag block_spread
            run 64 "$SPREAD" cyclic "${common[@]}" --tag cyclic_spread
            run 64 "$PACKED" block "${common[@]}" --tag block_packed
        done
    done
done
