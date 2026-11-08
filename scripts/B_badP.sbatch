#!/bin/bash
#SBATCH -J halo-B
#SBATCH -p eng-instruction
#SBATCH -N 4
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH -t 01:00:00
#SBATCH -o results/logs/%x-%j.out
# Experiment B: fixed global N, P = 48..64 with the most balanced grid for each P
# (MPI_Dims_create). Primes (53, 59, 61) are forced to 1x1xP.
cd "$SLURM_SUBMIT_DIR" && source scripts/common.sh

for N in 256 1024; do
    for P in $(seq 48 64); do
        for rep in $(seq "$REPS"); do
            run "$P" "$SPREAD" block --n "$N" --grid auto --steps "$(steps_for "$N")" \
                --mode nonblocking --pack manual --csv results/B.csv --tag B
        done
    done
done
