#!/bin/bash
#SBATCH -J halo-D
#SBATCH -p eng-instruction
#SBATCH -A 26fa-cs598pa-eng
#SBATCH --constraint=AE7713
#SBATCH -N 3
#SBATCH --ntasks-per-node=16
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH -t 01:00:00
#SBATCH -o results/logs/%x-%j.out
# Experiment D: blocking vs nonblocking vs overlap vs overlap+MPI_Test,
# on a good (4x4x3) and a bad (1x1x48) decomposition. Both transports: with rc
# only, every face above UCX_RNDV_THRESH (8 KB) uses the HCA rendezvous protocol,
# even within a node; with shm, intra-node faces use shared-memory rendezvous.
cd "$SLURM_SUBMIT_DIR" && source scripts/common.sh

for tls in shm rc; do
    set_tls "$tls"
    for N in 240 960; do
        for grid in 4x4x3 1x1x48; do
            for mode in blocking nonblocking overlap overlap_test; do
                for rep in $(seq "$REPS"); do
                    run 48 "$SPREAD" block --n "$N" --grid "$grid" --steps "$(steps_for "$N" 48)" \
                        --mode "$mode" --pack manual --csv results/D.csv --tag D
                done
            done
        done
    done
done
