#!/bin/bash
#SBATCH -J halo-ppmulti
#SBATCH -p eng-instruction
#SBATCH -A 26fa-cs598pa-eng
#SBATCH --constraint=AE7713
#SBATCH -N 2
#SBATCH --ntasks-per-node=16
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH -t 00:30:00
#SBATCH -o results/logs/%x-%j.out
# Bandwidth per pair with many pairs active at once (shm mode), to compare with
# the effective bandwidth of the A model. Rank r talks to rank r + P/2.
#   other node:   k = 1..16 pairs, cores 0,8,... on each node
#   other socket: 1 and 8 pairs on one node (cpu 8i <-> 8i+64)
#   same socket:  1 and 4 pairs on one node (cpu 8i <-> 8i+32)
cd "$SLURM_SUBMIT_DIR" && source scripts/common.sh
[ -x bin/pingpong_multi ] || { echo "bin/pingpong_multi missing: run make"; exit 1; }

pm() { local tag=$1 pat=$2; shift 2
    srun --whole "$@" bin/pingpong_multi --pattern "$pat" --csv results/pingpong_multi.csv --tag "$tag"; }
first() { seq 0 8 127 | head -n "$1" | paste -sd, -; }
for pat in pingpong exchange; do
    for k in 1 2 4 8 16; do
        pm other_node "$pat" -N 2 -n $((2 * k)) --ntasks-per-node=$k --cpu-bind=map_cpu:"$(first $k)"
    done
    pm other_socket "$pat" -N 1 -n 2 --cpu-bind=map_cpu:0,64
    pm other_socket "$pat" -N 1 -n 16 --cpu-bind=map_cpu:"$SPREAD"
    pm same_socket "$pat" -N 1 -n 2 --cpu-bind=map_cpu:0,32
    pm same_socket "$pat" -N 1 -n 8 --cpu-bind=map_cpu:"$(first 8)"
done
