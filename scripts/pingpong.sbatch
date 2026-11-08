#!/bin/bash
#SBATCH -J halo-pingpong
#SBATCH -p eng-instruction
#SBATCH -N 2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH -t 00:20:00
#SBATCH -o results/logs/%x-%j.out
# Alpha-beta measurement for each level of the hardware hierarchy.
cd "$SLURM_SUBMIT_DIR" && source scripts/common.sh

# pp TAG srun-placement-args...
pp() { local tag=$1; shift; srun -c 64 "$@" bin/pingpong --csv results/pingpong.csv --tag "$tag"; }
pp same_ccd     -N 1 -n 2 --cpu-bind=map_cpu:0,1
pp other_ccd    -N 1 -n 2 --cpu-bind=map_cpu:0,8
pp other_socket -N 1 -n 2 --cpu-bind=map_cpu:0,64
pp other_node   -N 2 -n 2 --ntasks-per-node=1 --cpu-bind=map_cpu:0
