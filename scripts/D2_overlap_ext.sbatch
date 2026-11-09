#!/bin/bash
#SBATCH -J halo-D2
#SBATCH -p eng-instruction
#SBATCH -A 26fa-cs598pa-eng
#SBATCH --constraint=AE7713
#SBATCH -N 3
#SBATCH --ntasks-per-node=16
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH -t 00:45:00
#SBATCH -o results/logs/%x-%j.out
# Experiment D2 (bin/halo3d_ext): completes D with N = 480 for the four original
# modes, and adds the peel mode (boundary first, then interior while the faces
# travel) with manual and indexed packing at N = 240, 480, 960.
cd "$SLURM_SUBMIT_DIR" && source scripts/common.sh
HALO_BIN=bin/halo3d_ext
[ -x $HALO_BIN ] || { echo "$HALO_BIN missing: run make"; exit 1; }

for tls in shm rc; do
    set_tls "$tls"
    for N in 240 480 960; do
        for grid in 4x4x3 1x1x48; do
            cases="peel:manual peel:indexed"
            [ "$N" = 480 ] && cases="blocking:manual nonblocking:manual overlap:manual overlap_test:manual $cases"
            for c in $cases; do
                for rep in $(seq "$REPS"); do
                    run 48 "$SPREAD" block --n "$N" --grid "$grid" --steps "$(steps_for "$N" 48)" \
                        --mode "${c%:*}" --pack "${c#*:}" --csv results/D2.csv --tag D2
                done
            done
        done
    done
done
