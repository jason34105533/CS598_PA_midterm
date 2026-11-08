#!/bin/bash
#SBATCH -J halo-C
#SBATCH -p eng-instruction
#SBATCH -A 26fa-cs598pa-eng
#SBATCH --constraint=AE7713
#SBATCH -N 3
#SBATCH --ntasks-per-node=16
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH -t 01:00:00
#SBATCH -o results/logs/%x-%j.out
# Experiment C: same grid, different rank placement.
#   block_spread : consecutive ranks on one node, one rank per CCD (baseline)
#   cyclic_spread: ranks dealt round-robin over nodes, one rank per CCD
#   block_packed : consecutive ranks on one node, 16 ranks packed on 2 CCDs
# Each under UCX_TLS=sm,self,rc (shm) and the site default UCX_TLS=rc (the tls
# column of C.csv tells them apart).
cd "$SLURM_SUBMIT_DIR" && source scripts/common.sh

for tls in shm rc; do
    set_tls "$tls"
    for N in 240 960; do
        for grid in 4x4x3 1x1x48 48x1x1; do
            for rep in $(seq "$REPS"); do
                common=(--n "$N" --grid "$grid" --steps "$(steps_for "$N" 48)" --mode nonblocking --pack manual --csv results/C.csv)
                run 48 "$SPREAD" block "${common[@]}" --tag block_spread
                run 48 "$SPREAD" cyclic "${common[@]}" --tag cyclic_spread
                run 48 "$PACKED" block "${common[@]}" --tag block_packed
            done
        done
    done
done
