#!/bin/bash
# Correctness sweep for bin/halo3d_ext: every mode x every packing (peel only with
# buffer packings), uneven subdomains and 1-cell-thick subdomains.
#   Local:   scripts/smoke_ext.sh
#   Cluster: LAUNCH="srun -N 2 -n" scripts/smoke_ext.sh   (inside an allocation)
LAUNCH=${LAUNCH:-"mpirun --oversubscribe -n"}
BIN=${BIN:-bin/halo3d_ext}

fail=0; count=0
check() { # check NRANKS args...
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

check 1 --n 12 --grid 1x1x1 --mode peel --pack indexed
check 7 --n 20 --grid auto --mode peel
check 8 --n 16 --grid 2x2x2 --mode peel --warmup 0
for grid in 2x2x2 8x1x1 1x8x1 1x1x8 2x4x1 1x2x4; do
    for N in 9 23; do
        for mode in blocking nonblocking overlap overlap_test peel; do
            for pack in manual datatype indexed indexed_dt; do
                [[ $mode == peel && ( $pack == datatype || $pack == indexed_dt ) ]] && continue
                check 8 --n "$N" --grid "$grid" --mode "$mode" --pack "$pack"
            done
        done
    done
done

echo "smoke_ext: $((count - fail))/$count passed"
[ "$fail" -eq 0 ]
