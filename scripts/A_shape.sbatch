#!/bin/bash
#SBATCH -J halo-A
#SBATCH -p eng-instruction
#SBATCH -A 26fa-cs598pa-eng
#SBATCH --constraint=AE7713
#SBATCH -N 3
#SBATCH --ntasks-per-node=16
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH -t 01:30:00
#SBATCH -o results/logs/%x-%j.out
# Experiment A: fixed P=48, every factorization Px x Py x Pz (orientation included),
# small and large problem.
cd "$SLURM_SUBMIT_DIR" && source scripts/common.sh

for N in 240 960; do
    for grid in $(factor3 48); do
        for rep in $(seq "$REPS"); do
            run 48 "$SPREAD" block --n "$N" --grid "$grid" --steps "$(steps_for "$N" 48)" \
                --mode nonblocking --pack manual --csv results/A.csv --tag A
        done
    done
done
