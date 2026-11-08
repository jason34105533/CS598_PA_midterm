#!/bin/bash
#SBATCH -J halo-S2
#SBATCH -p eng-instruction
#SBATCH -A 26fa-cs598pa-eng
#SBATCH --constraint=AE7713
#SBATCH -N 3
#SBATCH --ntasks-per-node=16
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH -t 00:30:00
#SBATCH -o results/logs/%x-%j.out
# Experiment S2: S with P <= 16 placed alternately on the two sockets (ALT map)
# instead of filling socket 0 first (SPREAD). Same N and steps as S; the P = 1
# baseline comes from S. Tags strong_alt / weak240_alt in results/S.csv.
cd "$SLURM_SUBMIT_DIR" && source scripts/common.sh
echo "ALT=$ALT"

PS="2 3 4 6 8 12 16"
for N in 480 960; do
    for P in $PS; do
        for rep in $(seq "$REPS"); do
            run "$P" "$ALT" block --n "$N" --grid auto --steps "$(steps_for "$N" "$P")" \
                --mode nonblocking --pack manual --csv results/S.csv --tag strong_alt
        done
    done
done
for P in $PS; do
    N=$(awk -v p="$P" 'BEGIN { printf "%d", 240 * p ^ (1 / 3) + 0.5 }')
    for rep in $(seq "$REPS"); do
        run "$P" "$ALT" block --n "$N" --grid auto --steps "$(steps_for "$N" "$P")" \
            --mode nonblocking --pack manual --csv results/S.csv --tag weak240_alt
    done
done
