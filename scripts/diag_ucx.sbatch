#!/bin/bash
#SBATCH -J halo-diag
#SBATCH -p eng-instruction
#SBATCH -A 26fa-cs598pa-eng
#SBATCH --constraint=AE7713
#SBATCH -N 2
#SBATCH --ntasks-per-node=2
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH -t 00:15:00
#SBATCH -o results/logs/%x-%j.out
# Does intra-node traffic above 8 KB use shared memory?
#   shm      : UCX_TLS=sm,self,rc (rendezvous above 8 KB everywhere)
#   shm_eager: same, plus UCX_RNDV_THRESH=intra:inf,inter:8K (intra-node always eager)
cd "$SLURM_SUBMIT_DIR" && source scripts/common.sh

echo "== compute node: ptrace_scope $(srun -N 1 -n 1 cat /proc/sys/kernel/yama/ptrace_scope)"
srun -N 1 -n 1 ucx_info -d | grep -E 'Transport:' | sort | uniq -c

pp() { local tag=$1; shift; srun --whole "$@" bin/pingpong --max-bytes 4194304 --csv results/diag_ucx.csv --tag "$tag"; }
for cfg in shm shm_eager; do
    set_tls shm
    [ "$cfg" = shm_eager ] && export UCX_RNDV_THRESH=intra:inf,inter:8K
    echo "== $cfg: UCX_TLS=$UCX_TLS UCX_RNDV_THRESH=$UCX_RNDV_THRESH"
    UCX_PROTO_INFO=y pp "${cfg}_same_ccd" -N 1 -n 2 --cpu-bind=map_cpu:0,1 2>&1 | grep -E 'ucp_ep|rendezvous|eager|[0-9]+ +[0-9.]+ +[0-9.]+' | head -60
    pp "${cfg}_other_socket" -N 1 -n 2 --cpu-bind=map_cpu:0,64
    pp "${cfg}_other_node"   -N 2 -n 2 --ntasks-per-node=1 --cpu-bind=map_cpu:0
done
