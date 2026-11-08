#!/bin/bash
#SBATCH -J halo-D
#SBATCH -p eng-instruction
#SBATCH -N 4
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH -t 01:00:00
#SBATCH -o results/logs/%x-%j.out
# Experiment D: blocking vs nonblocking vs overlap vs overlap+MPI_Test,
# on a good (4x4x4) and a bad (1x1x64) decomposition.
cd "$SLURM_SUBMIT_DIR" && source scripts/common.sh

for N in 256 1024; do
    for grid in 4x4x4 1x1x64; do
        for mode in blocking nonblocking overlap overlap_test; do
            for rep in $(seq "$REPS"); do
                run 64 "$SPREAD" block --n "$N" --grid "$grid" --steps "$(steps_for "$N")" \
                    --mode "$mode" --pack manual --csv results/D.csv --tag D
            done
        done
    done
done
