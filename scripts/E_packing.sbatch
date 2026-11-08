#!/bin/bash
#SBATCH -J halo-E
#SBATCH -p eng-instruction
#SBATCH -N 4
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH -t 01:00:00
#SBATCH -o results/logs/%x-%j.out
# Experiment E: manual packing vs MPI derived datatypes. The 1-D grids isolate
# one face direction each (x: strided, y: rows, z: rows of full planes).
cd "$SLURM_SUBMIT_DIR" && source scripts/common.sh

for N in 256 1024; do
    for grid in 64x1x1 1x64x1 1x1x64 4x4x4; do
        for pack in manual datatype; do
            for rep in $(seq "$REPS"); do
                run 64 "$SPREAD" block --n "$N" --grid "$grid" --steps "$(steps_for "$N")" \
                    --mode nonblocking --pack "$pack" --csv results/E.csv --tag E
            done
        done
    done
done
