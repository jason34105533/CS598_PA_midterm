#!/bin/bash
#SBATCH -J halo-B
#SBATCH -p eng-instruction
#SBATCH -A 26fa-cs598pa-eng
#SBATCH --constraint=AE7713
#SBATCH -N 3
#SBATCH --ntasks-per-node=16
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH -t 01:00:00
#SBATCH -o results/logs/%x-%j.out
# Experiment B: fixed global N, P = 36..48 with the most balanced grid for each P
# (MPI_Dims_create). P < 48 spreads ranks evenly over 3 nodes. Primes (37, 41, 43, 47) are forced to 1x1xP.
cd "$SLURM_SUBMIT_DIR" && source scripts/common.sh

for N in 240 960; do
    for P in $(seq 36 48); do
        for rep in $(seq "$REPS"); do
            run "$P" "$SPREAD" block --n "$N" --grid auto --steps "$(steps_for "$N" "$P")" \
                --mode nonblocking --pack manual --csv results/B.csv --tag B
        done
    done
done
