#!/bin/bash
# Correctness sweep over grid shapes x modes x packing, with N not divisible
# by the grid (uneven subdomains) and subdomains only 1 cell thick.
#   Local:   scripts/smoke.sh
#   Cluster: LAUNCH="srun -N 2 -n" scripts/smoke.sh   (inside an allocation)
LAUNCH=${LAUNCH:-"mpirun --oversubscribe -n"}
BIN=${BIN:-bin/halo3d}

fail=0; count=0
check() { # check NRANKS halo3d-args...
    local np=$1; shift
    local out
    out=$($LAUNCH "$np" "$BIN" --steps 7 --warmup 2 "$@" 2>&1)
    count=$((count + 1))
    if ! grep -q PASS <<<"$out"; then
        fail=$((fail + 1))
        echo "FAIL: np=$np $*"
        echo "$out" | tail -5
    fi
}

check 1 --n 12 --grid 1x1x1
check 7 --n 20 --grid auto
for grid in 2x2x2 8x1x1 1x8x1 1x1x8 2x4x1 1x2x4; do
    for N in 9 23; do
        for mode in blocking nonblocking overlap overlap_test; do
            for pack in manual datatype; do
                check 8 --n "$N" --grid "$grid" --mode "$mode" --pack "$pack"
            done
        done
    done
done

echo "smoke: $((count - fail))/$count passed"
[ "$fail" -eq 0 ]
