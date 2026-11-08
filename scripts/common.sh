# Sourced by every job script (run from the repo root).

# Edit the module names to match `module avail openmpi` on Campus Cluster
# (HW3 used gcc 13.3 + Open MPI 5.0.1 over UCX).
module load gcc openmpi 2>/dev/null
command -v mpicc >/dev/null || { echo "mpicc not found: fix module names in scripts/common.sh"; exit 1; }
[ -x bin/halo3d ] || { echo "bin/halo3d missing: run 'make ARCH=-march=znver3' on the login node"; exit 1; }

mkdir -p results/logs
REPS=${REPS:-3}

# 2x EPYC 7713: 128 cores, 16 CCDs of 8 consecutive cores (check with lscpu, see smoke.sbatch).
SPREAD=$(seq 0 8 127 | paste -sd, -)  # one rank per CCD: 0,8,...,120
PACKED=$(seq 0 15 | paste -sd, -)      # 16 ranks on 2 CCDs

# Steps per run, chosen so one run takes about 1-2 s.
steps_for() { case $1 in 256) echo 1000 ;; 512) echo 300 ;; *) echo 100 ;; esac; }

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
# -c 8 makes each step own whole nodes (16 tasks x 8 CPUs), so map_cpu can pick
# any core; srun does not inherit --cpus-per-task from sbatch on recent Slurm.
run() {
    local n=$1 map=$2 dist=$3
    shift 3
    srun -n "$n" -c 8 --distribution="$dist" --cpu-bind=map_cpu:"$map" bin/halo3d "$@"
}

echo "job $SLURM_JOB_ID on $SLURM_JOB_NODELIST  ($(date))"
env | grep -E '^(UCX|OMPI_MCA)_' | sort
