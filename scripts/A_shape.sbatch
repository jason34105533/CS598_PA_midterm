#!/bin/bash
#SBATCH -J halo-A
#SBATCH -p eng-instruction
#SBATCH -N 4
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH -t 01:30:00
#SBATCH -o results/logs/%x-%j.out
# Experiment A: fixed P=64, every factorization Px x Py x Pz (orientation included),
# small and large problem.
cd "$SLURM_SUBMIT_DIR" && source scripts/common.sh

for N in 256 1024; do
    for grid in $(factor3 64); do
        for rep in $(seq "$REPS"); do
            run 64 "$SPREAD" block --n "$N" --grid "$grid" --steps "$(steps_for "$N")" \
                --mode nonblocking --pack manual --csv results/A.csv --tag A
        done
    done
done
