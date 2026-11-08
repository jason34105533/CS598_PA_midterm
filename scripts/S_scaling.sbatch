#!/bin/bash
#SBATCH -J halo-S
#SBATCH -p eng-instruction
#SBATCH -A 26fa-cs598pa-eng
#SBATCH --constraint=AE7713
#SBATCH -N 3
#SBATCH --ntasks-per-node=16
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH -t 01:30:00
#SBATCH -o results/logs/%x-%j.out
# Experiment S: strong and weak scaling with the auto grid, one rank per CCD,
# ceil(P/16) nodes. Gives parallel efficiency vs n/P and n_0.8.
#   strong: fixed N, P = 1..48
#   weak:   N = round(L * P^(1/3)), so the local volume stays near L^3
#           (exact only when P is a cube; cells_min/cells_max record the spread)
cd "$SLURM_SUBMIT_DIR" && source scripts/common.sh

PS="1 2 3 4 6 8 12 16 24 32 48"
for N in 240 480 960; do
    for P in $PS; do
        for rep in $(seq "$REPS"); do
            run "$P" "$SPREAD" block --n "$N" --grid auto --steps "$(steps_for "$N" "$P")" \
                --mode nonblocking --pack manual --csv results/S.csv --tag strong
        done
    done
done
for L in 60 240; do
    for P in $PS; do
        N=$(awk -v l="$L" -v p="$P" 'BEGIN { printf "%d", l * p ^ (1 / 3) + 0.5 }')
        for rep in $(seq "$REPS"); do
            run "$P" "$SPREAD" block --n "$N" --grid auto --steps "$(steps_for "$N" "$P")" \
                --mode nonblocking --pack manual --csv results/S.csv --tag "weak$L"
        done
    done
done
