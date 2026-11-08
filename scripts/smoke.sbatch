#!/bin/bash
#SBATCH -J halo-smoke
#SBATCH -p eng-instruction
#SBATCH -A 26fa-cs598pa-eng
#SBATCH --constraint=AE7713
#SBATCH -N 3
#SBATCH --ntasks-per-node=16
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH -t 00:20:00
#SBATCH -o results/logs/%x-%j.out
# Correctness on the cluster + hardware/placement sanity checks.
cd "$SLURM_SUBMIT_DIR" && source scripts/common.sh

echo "== CPU -> core/socket/NUMA/L3 (CCD = cores sharing one L3)"
srun -N 1 -n 1 lscpu | grep -E 'Model name|Socket|NUMA node\(s\)|L3'
srun -N 1 -n 1 lscpu -e=CPU,CORE,SOCKET,NODE,CACHE

echo "== correctness sweep across 2 nodes"
LAUNCH="srun -N 2 -n" scripts/smoke.sh || exit 1

echo "== placement check: P=48, one rank per CCD"
run 48 "$SPREAD" block --n 240 --grid 4x4x3 --steps 20 --verbose
