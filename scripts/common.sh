# Sourced by every job script (run from the repo root).

# gcc/13.3.0 + openmpi/5.0.1-gcc-13.3.0 (UCX) on Campus Cluster.
module load gcc openmpi 2>/dev/null
command -v mpicc >/dev/null || { echo "mpicc not found: fix module names in scripts/common.sh"; exit 1; }
[ -x bin/halo3d ] || { echo "bin/halo3d missing: run 'make ARCH=-march=znver3' on the login node"; exit 1; }

# Without this, srun starts every rank as an MPI singleton (size 1).
export SLURM_MPI_TYPE=pmix

# UCX transports. The site default is UCX_TLS=rc, which sends even intra-node
# messages through the InfiniBand HCA (loopback): about 1.3 us and 8 GB/s within
# one CCD. "shm" adds the shared-memory transports, so intra-node messages stay
# in memory and only inter-node ones use rc. All experiments use shm unless
# they loop over both (pingpong, C).
#   set_tls shm|rc
set_tls() {
    case $1 in
        shm) export UCX_TLS=sm,self,rc ;;
        rc)  export UCX_TLS=rc ;;
        *)   echo "set_tls: unknown mode $1"; exit 1 ;;
    esac
}
set_tls "${TLS_MODE:-shm}"

mkdir -p results/logs
REPS=${REPS:-3}

# 2x EPYC 7713: 128 cores, 16 CCDs of 8 consecutive cores (check with lscpu, see smoke.sbatch).
SPREAD=$(seq 0 8 127 | paste -sd, -)  # one rank per CCD: 0,8,...,120
PACKED=$(seq 0 15 | paste -sd, -)      # 16 ranks on 2 CCDs

# steps_for N P: about 1e9 cell updates per rank (roughly 1 s), clamped to [20, 1000].
steps_for() {
    local s=$((1000000000 * $2 / ($1 * $1 * $1)))
    ((s < 20)) && s=20
    ((s > 1000)) && s=1000
    echo "$s"
}

# All ordered factorizations P = Px*Py*Pz, printed as PxxPyxPz.
factor3() {
    local P=$1 px py r
    for ((px = 1; px <= P; px++)); do
        ((P % px)) && continue
        r=$((P / px))
        for ((py = 1; py <= r; py++)); do
            ((r % py)) && continue
            echo "${px}x${py}x$((r / py))"
        done
    done
}

# run NTASKS CPU_MAP DISTRIBUTION halo3d-args...
# Uses ceil(NTASKS/16) nodes. --whole gives the step every CPU of its nodes, so
# map_cpu can pick any core even when NTASKS < 16 per node. --whole is ignored
# when -c is set, which is why the job headers do not use --cpus-per-task.
# With NTASKS not a multiple of 16, Slurm spreads tasks evenly (37 -> 13/12/12).
run() {
    local n=$1 map=$2 dist=$3
    shift 3
    srun -n "$n" -N $(((n + 15) / 16)) --ntasks-per-node=16 --whole \
        --distribution="$dist" --cpu-bind=map_cpu:"$map" bin/halo3d "$@"
}

echo "job $SLURM_JOB_ID on $SLURM_JOB_NODELIST  ($(date))"
env | grep -E '^(UCX|OMPI_MCA)_' | sort
