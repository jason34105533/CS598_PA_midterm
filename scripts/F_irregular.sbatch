#!/bin/bash
#SBATCH -J halo-F
#SBATCH -p eng-instruction
#SBATCH -A 26fa-cs598pa-eng
#SBATCH --constraint=AE7713
#SBATCH -N 3
#SBATCH --ntasks-per-node=16
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH -t 00:30:00
#SBATCH -o results/logs/%x-%j.out
# Experiment F: layered irregular decomposition (bin/halo3d_irr) for prime P,
# against 1x1xP and against using fewer ranks with a regular grid (bin/halo3d,
# auto grid, same settings as B). 48 = 4x4,4x4,4x4 checks that halo3d_irr matches
# halo3d on a regular grid.
cd "$SLURM_SUBMIT_DIR" && source scripts/common.sh
[ -x bin/halo3d_irr ] || { echo "bin/halo3d_irr missing"; exit 1; }

declare -A LAYERS=([37]=4x4,3x4,3x3 [41]=4x4,4x4,3x3 [43]=4x4,3x5,3x4 [47]=4x4,4x4,3x5 [48]=4x4,4x4,4x4)
irr() { # irr NTASKS args...
    local n=$1; shift
    srun -n "$n" -N $(((n + 15) / 16)) --ntasks-per-node=16 --whole \
        --distribution=block --cpu-bind=map_cpu:"$SPREAD" bin/halo3d_irr "$@"
}

for N in 240 960; do
    for P in 37 41 43 47 48; do
        for rep in $(seq "$REPS"); do
            irr "$P" --n "$N" --layers "${LAYERS[$P]}" --steps "$(steps_for "$N" "$P")" \
                --csv results/F.csv --tag irr
        done
    done
    for P in 36 37 40 41 42 43 45 47 48; do
        for rep in $(seq "$REPS"); do
            run "$P" "$SPREAD" block --n "$N" --grid auto --steps "$(steps_for "$N" "$P")" \
                --mode nonblocking --pack manual --csv results/F_ref.csv --tag ref
        done
    done
done
