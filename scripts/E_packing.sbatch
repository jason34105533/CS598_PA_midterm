#!/bin/bash
#SBATCH -J halo-E
#SBATCH -p eng-instruction
#SBATCH -A 26fa-cs598pa-eng
#SBATCH --constraint=AE7713
#SBATCH -N 3
#SBATCH --ntasks-per-node=16
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH -t 01:00:00
#SBATCH -o results/logs/%x-%j.out
# Experiment E: manual packing vs MPI derived datatypes. The 1-D grids isolate
# one face direction each (x: strided, y: rows, z: rows of full planes).
cd "$SLURM_SUBMIT_DIR" && source scripts/common.sh

for N in 240 960; do
    for grid in 48x1x1 1x48x1 1x1x48 4x4x3; do
        for pack in manual datatype; do
            for rep in $(seq "$REPS"); do
                run 48 "$SPREAD" block --n "$N" --grid "$grid" --steps "$(steps_for "$N" 48)" \
                    --mode nonblocking --pack "$pack" --csv results/E.csv --tag E
            done
        done
    done
done
