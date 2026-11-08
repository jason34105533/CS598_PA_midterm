# CS598 PA Midterm Project Plan: 3D Halo Exchange and Domain Decomposition

Team: Matei Budiu, Wei-Sheng Lin. Due: Nov 11, 2026 (presentation + short written document).

## 1. Storyline

The domain decomposition decides four things about communication:

- how much data moves (halo surface area),
- how many messages are sent (number of neighbors),
- whether that data is contiguous in memory (packing cost),
- which level of the hardware it crosses (core, CCD, node).

The other techniques (overlap, packing, rank placement) either hide or amplify what the decomposition causes. Every experiment below is organized around this idea.

## 2. Model problem

3D heat diffusion with explicit time stepping on the unit cube, with Dirichlet zero boundaries:

```
u_new[i,j,k] = u[i,j,k] + r * ( u[i-1,j,k] + u[i+1,j,k]
                              + u[i,j-1,k] + u[i,j+1,k]
                              + u[i,j,k-1] + u[i,j,k+1] - 6*u[i,j,k] )
```

- Each cell reads its 6 face neighbors (7-point stencil), so each rank needs one ghost layer per face. Only faces are exchanged; edges and corners are not needed.
- This is the core pattern of structured-mesh CFD/diffusion codes (OpenFOAM-style "update interfaces, compute internal, compute boundary"), with all other physics removed.
- A fixed number of time steps is used. There is no convergence study.
- `r = 1/8` (stable, since `r <= 1/6` is required).

**Correctness check.** With the initial condition `u0 = sin(pi x) sin(pi y) sin(pi z)`, the discrete solution after `k` steps is exactly `g^k * u0`, where `g = 1 - 12 r sin^2(pi h / 2)`. The max error against this should be at machine precision for every decomposition and every communication mode. A checksum against a P=1 run serves as a second check.

## 3. What "bad decomposition" means

Three separate cases:

1. **Bad shape at fixed P.** For example, 480 = 1x1x480 or 6x8x10. The compute per rank is the same, but the halo surface differs a lot. Near-cubic subdomains minimize the surface.
2. **Bad P.** If P is prime or has a large prime factor, the only option may be 1x1xP (for example P=61). Using slightly fewer ranks (60 = 3x4x5) costs about 1.6% more work per rank but roughly 5x less halo.
3. **N not divisible by the grid.** This causes uneven subdomains (ceil/floor) or padding. A stretch idea is an irregular decomposition that beats padding or the 1x1xP fallback.

Halo cells per interior rank, for global `N^3` on a `Px x Py x Pz` grid:

```
halo = 2 N^2 * ( [Px>1]/(Py*Pz) + [Py>1]/(Px*Pz) + [Pz>1]/(Px*Py) )
```

| Grid (P=48) | Neighbors | Halo per rank (units of N^2) |
|---|---|---|
| 1x1x48 | 2 | 2.0 |
| 2x4x6 | 6 | 0.50 |
| 4x4x3 | 6 | 0.46 |

The slab moves about 4.4x more data (roughly P^(2/3)/3) but sends 3x fewer messages, so it may win when messages are small and latency dominates. Compare against the alpha-beta model `T = m*alpha + bytes*beta`.

**Shape also changes memory behavior, not only the surface area.** Data is stored with x as the unit-stride direction, so:

- Splitting along x shortens inner loops, which hurts vectorization and prefetch. It also makes x-faces strided, which makes packing expensive.
- Splitting only along z (1x1xP) makes z-faces contiguous, so no packing is needed. However, each xy-plane is large, so the 3 planes the stencil needs may not fit in L2.
- So 1x1x48 and 48x1x1 have the same surface, but performance may differ. Orientation is an experiment in its own right.

**Hierarchy: where does the extra surface land?** Take 3 nodes with 16 ranks per node, ranks filled into nodes in order, and x varying fastest in rank numbering:

- 1x1x48: each node owns 16 consecutive slabs, so the off-node faces are two N^2 planes.
- 4x4x3: each node owns one 4x4x1 layer, so the off-node faces are also two N^2 planes.

Inter-node traffic is identical, because the node-level decomposition is 1x1x3 in both cases. The extra surface of the slab is all intra-node. Whether intra-node is actually cheaper depends on the MPI transport (Section 8 and the measurement log): with the site default `UCX_TLS=rc`, intra-node messages also go through the InfiniBand HCA. This is why 4x4x3 is the "good" grid in C, D and E rather than 3x4x4 (same surface, but a z-layer of 12 ranks does not line up with a 16-rank node). With cyclic rank distribution, the same 1x1x48 sends every face off-node, which is the worst case. So decomposition quality has to be judged at each level: node, CCD, and core. The CCD level comes in naturally here.

## 4. Code design

One C + MPI program. Every variant is a command-line option.

- Process grid `Px Py Pz` is given explicitly (or chosen by `MPI_Dims_create`). Uneven splits are allowed.
- Local array with one ghost layer. Global size `N` (or local size) is configurable.
- Exchange modes:
  - `blocking`: `MPI_Sendrecv` per direction, then compute everything.
  - `nonblocking`: `Isend/Irecv` on all faces, `Waitall`, then compute everything.
  - `overlap`: `Isend/Irecv`, compute the interior, `Waitall`, compute the boundary layer.
  - `overlap_test`: same as `overlap`, with periodic `MPI_Test` calls during the interior compute (forces progress of rendezvous messages).
- Packing: `manual` (pack into contiguous buffers) or `datatype` (`MPI_Type_create_subarray` / `MPI_Type_vector`).
- Per-phase timers: pack, post, interior compute, wait, boundary compute, unpack, total. Report max and average over ranks.
- Correctness check from Section 2, run for every mode.
- Rank placement is controlled only through `srun` (block/cyclic distribution, `--cpu-bind=map_cpu`). The code does not change.
- Helper: a small ping-pong test to measure alpha and beta between cores in the same CCD, across CCDs, and across nodes.

## 5. Experiments (in priority order)

| ID | Experiment | Question |
|---|---|---|
| A | Fixed P, sweep all factorizations including orientation (1x1xP vs Px1x1). Two problem sizes: large and small local volume. | How much does shape matter? Does the surface/alpha-beta model predict it? Is there a latency/bandwidth crossover? |
| S | Strong scaling (fixed N, P = 1..48) and weak scaling (local volume about L^3), auto grid. | Parallel efficiency vs n/P; what n/P gives 80% efficiency (n_0.8)? Is the loss compute (cache) or communication? |
| B | Prime P vs slightly fewer ranks with a good factorization (for example 47 vs 48), fixed global N. | Can using fewer cores be faster? |
| C | Same grid, different placement: block vs cyclic across nodes; inside a node, one rank per CCD vs packed into few CCDs. Each under `shm` and `rc` transport. | Does it matter where the extra surface lands (node / CCD / core)? Does the answer change when intra-node messages go through the NIC? |
| D | blocking vs overlap vs overlap_test, on the best and worst decomposition from A, under `shm` and `rc`. | How much communication can overlap hide? Does it really overlap (compare Waitall time with and without `MPI_Test`)? How does intra-node rendezvous (HCA vs shared memory) change this? |
| E | manual packing vs derived datatypes, broken down by face direction. | How expensive are non-contiguous faces? Does the answer depend on the decomposition? |
| F (stretch) | Irregular decomposition for prime P. | Can it beat 1x1xP and beat dropping ranks? |

## 5a. Predictions (written before the main runs)

Common setup: 3 nodes x 16 ranks (one per CCD), `shm` transport unless stated, 3 repetitions, median reported. Main metric `t_step` (slowest rank), plus the phase breakdown and halo bytes per hierarchy level. These are expectations to check against the data, not results.

**A thread through all experiments: the two-point alpha-beta model probably underestimates communication.** In the smoke run (N=240, 4x4x3, `rc`), `wait` was 0.64 ms of a 0.87 ms step, while pingpong alpha/beta predicts about 30 us for the 211 KB each rank sends. A likely cause is that the 16 ranks of a node share one HCA (under `rc` even intra-node traffic uses it), so per-node injection bandwidth, not the per-pair beta, sets the cost. The model will probably need a per-node bandwidth term: `T_comm ~ max(sum over a rank's messages of (alpha + m*beta), node_bytes / B_node)`.

**S (scaling).**
- Strong: N = 240, 480, 960; P = 1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48. Weak: local volume about 60^3 and 240^3.
- Plots: time vs P (log-log), efficiency vs n/P (semi-log x, one curve per N), n_0.8.
- Expect a kink between P = 16 (one node) and P = 24 (first inter-node traffic).
- N = 480 at P = 48 has two arrays of about 37 MB per rank, right at the 32 MB L3 of a CCD, so superlinear compute speedup is possible there. That is the "cache or communication" question in concrete form.
- N = 240 at large P should be communication-bound with low efficiency.

**A (shape, the core experiment).** All 45 ordered factorizations of 48, N = 240, 480 and 960 (480 puts two arrays per rank at about 37 MB, at the edge of the 32 MB L3). Group them as 1D (slab), 2D (pencil), 3D (block) decompositions.
- N = 240 (communication-bound): 1x1x48 sends 922 KB per rank in 2 messages, 4x4x3 sends 211 KB in 6. With alpha around 1 us, the latency term is small, so the slab should be clearly slower. A latency/bandwidth crossover would need a much smaller N.
- N = 960 (compute-bound): compute is about 20 ms per step per rank, communication 1-2 ms, so shape matters mostly through memory behavior. 48x1x1 leaves 20 cells in the unit-stride direction (short inner loops, poor vectorization and prefetch) and strided x-faces. 1x1x48 and 48x1x1 have the same surface but may differ a lot (orientation effect).
- Analysis: measured `t_step` against halo bytes and message count, with the model `T = T_comp + T_comm` (including the node-bandwidth term) on top.

**B (bad P).** N fixed (240, 960), P = 36..48 with `MPI_Dims_create`. 37, 41, 43, 47 are forced to 1x1xP.
- N = 240: P = 47 (1x1x47) should be much slower than P = 48 (3x4x4), and possibly slower than P = 40 (2x4x5) or 45 (3x3x5), which would be the "fewer cores is faster" case.
- N = 960: compute dominates, so P = 47 should be only about 2% slower than 48.
- Caveat: for P not a multiple of 16, Slurm spreads ranks evenly (37 -> 13/12/12), so node placement differs between P values. Use `bytes_offnode` to separate that effect.

**C (placement x transport).** Grids 4x4x3, 1x1x48, 48x1x1; block_spread (baseline), cyclic_spread (every slab face crosses nodes), block_packed (16 ranks of a node on 2 CCDs); each under `shm` and `rc`.
- `shm`: the block vs cyclic gap should be largest, because intra-node traffic is cheap and cyclic pushes all of it onto the network.
- `rc`: intra-node traffic also goes through the HCA, so the block vs cyclic gap should shrink. The hierarchy argument of Section 3 only holds if the levels really differ in cost, and the transport decides that.
- packed: communication is closer, but 8 ranks share one L3 and the CCD's memory bandwidth, so compute at N = 960 should slow down. This shows directly why undersubscription is used.

**D (overlap x transport).** Modes `blocking`, `nonblocking`, `overlap`, `overlap_test`; grids 4x4x3 and 1x1x48; `shm` and `rc`.
- Mode names: all four exchange the same halos and differ only in where compute goes. `nonblocking` posts Isend/Irecv and immediately calls Waitall, then computes, so it does not overlap anything. It is the control. `overlap` puts the interior compute between the post and the Waitall.
- Every face is above `UCX_RNDV_THRESH` (8 KB), so it uses the rendezvous protocol (RTS -> RTR -> RDMA put -> FIN). Each handshake step is processed only while one side is inside an MPI call, and UCX has no progress thread by default. During the interior compute neither side calls MPI, so the handshake is expected to stall until Waitall. Prediction: `overlap` has about the same `wait` as `nonblocking` (no real overlap), and `overlap_test` (periodic `MPI_Testall`) shows a clearly smaller `wait`. Once the RDMA starts, the HCA finishes it without the CPU, so a few Test calls should be enough. Eager messages (below 8 KB) would overlap without Test.
- The hideable amount is bounded by the interior compute: at N = 240, 0.18 ms of compute against 0.64 ms of communication; at N = 960, compute is much larger than communication.
- The cost of overlap is the separate boundary pass (`comp_bd`) and worse cache use. Check whether it eats the gain.
- `shm` vs `rc`: intra-node rendezvous uses shared memory in one case and HCA RDMA in the other, which may depend on progress differently.

**E (packing).** Grids 48x1x1 (x-faces only, fully strided, worst case), 1x48x1 (y-faces, one contiguous row block per z-plane), 1x1x48 (z-faces, contiguous), 4x4x3 (mixed). Manual vs datatype.
- z-faces: the datatype path might avoid the copy (zero-copy) and beat manual packing.
- x-faces: expensive either way, and the datatype path may be slower than manual packing.
- Plot pack/unpack time against wait time, by face direction.

**F (stretch).** An irregular decomposition for prime P (for example, 47 as 48 blocks with one rank taking two), compared with 1x1x47 and with dropping one rank.

**Storyline.** A shows that shape sets the communication volume and the memory behavior. B covers the case where P limits the shape. C shows that the same volume costs differently depending on which level it lands on, and that the transport sets those costs. D and E ask whether communication can be hidden or made cheaper. S provides the scaling and efficiency numbers.

## 6. Methodology

- Baseline: multiple nodes, undersubscribed, one rank per CCD (16 ranks per node), as the instructor suggested.
- Fixed number of time steps, a few warmup steps, `MPI_Barrier` before timing.
- Each configuration is repeated at least 3 times; the median is reported.
- Platform: UIUC Campus Cluster, eng-instruction partition, 2x AMD EPYC 7713 per node (128 cores, no SMT, 16 CCDs of 8 cores, 32 MB L3 per CCD), HDR InfiniBand. gcc 13.3.0, Open MPI 5.0.1 over UCX. Record UCX settings (for example the rendezvous threshold, which matters for D).

## 7. Schedule and risk

- The real bottleneck is the cluster queue, not coding. Have the code correct early and batch many configurations into a few jobs.
- Nov 7: code written; correctness sweep passes locally.
- Nov 8: cluster smoke test, then submit all experiments at once (`scripts/submit_all.sh`).
- Nov 9: plots and analysis; rerun anything suspicious.
- Nov 10: written document and slides.

## 8. Decisions

- Partition: eng-instruction, account `26fa-cs598pa-eng`, `--constraint=AE7713`. The partition has 5 nodes, but only 3 (ccc0391-0393) are EPYC 7713; the other 2 are 64-core Intel Xeon 8358 nodes, so every job uses those 3 nodes (12 h limit, exclusive). `secondary` has more 7713 nodes (`--constraint=AE7713_100g_512G_NoGPU`) but a 4 h limit, shared nodes and a long queue.
- Baseline: 3 nodes x 16 ranks (one per CCD) = P = 48 (multiple nodes, undersubscribed). Experiment B uses P = 36..48, which includes the primes 37, 41, 43 and 47. For P < 48, Slurm spreads ranks evenly over the nodes (37 -> 13/12/12).
- Problem sizes: N = 240 (local 60x60x80 at 4x4x3, about 2.3 MB per array, fits in L3, so communication matters more) and N = 960 (local 240x240x320, about 147 MB per array, memory-bound compute dominates). 240 and 960 are divisible by every factor of 48, so the A/C/D/E grids have even subdomains; uneven splits show up only in B and S.
- Steps per run: about 1e9 cell updates per rank, clamped to 20..1000 (`steps_for` in `scripts/common.sh`).
- MPI transport: the site environment sets `UCX_TLS=rc` (plus `OMPI_MCA_pml=ucx`, `UCX_RNDV_THRESH=8192`, `UCX_RNDV_SCHEME=put_zcopy`). With rc only, every message, including one between two cores of the same CCD, goes through the HCA. All experiments use `shm` (`UCX_TLS=sm,self,rc` plus `UCX_RNDV_THRESH=intra:inf,inter:8K`: every intra-node message is an eager copy through shared memory, inter-node messages use rc). The threshold is needed because CMA is disabled on the nodes (`ptrace_scope=3`), so without it intra-node rendezvous still goes through the HCA. pingpong, C and D run under both, as an extra dimension (`set_tls shm|rc` in `scripts/common.sh`; the `tls` CSV column records it).
- Launch: every `srun` needs PMIx (`SLURM_MPI_TYPE=pmix`, set in `scripts/common.sh`); without it each rank starts as a singleton.

## 8a. Measurement log

**2026-11-08, job 12258647 (smoke, `UCX_TLS=rc`).** 98/98 correctness checks passed across 2 nodes. Placement check at P=48, 4x4x3, N=240: ranks 0-15 on ccc0391 at cpus 0,8,...,120 (one per CCD), and each node holds one z-layer. Off-node halo traffic is 1.84 MB per step, which matches two 240^2 planes in each direction. 20 steps: 0.87 ms per step, of which `wait` is 0.64 ms and compute 0.18 ms.

**2026-11-08, job 12258648 (pingpong v0, `UCX_TLS=rc`).** Mean of one trial per size; mid sizes are noisy, so this version was replaced by min/median over 7 trials. Raw data: `results/pingpong_v0_rc.csv`.

| Pair | 8 B (us) | 8 KB (us) | 1 MB (GB/s) | 16 MB (GB/s) |
|---|---|---|---|---|
| same CCD (cpu 0,1) | 1.37 | 10.5 | 7.7 | 8.0 |
| other CCD (cpu 0,8) | 1.35 | 9.8 | 7.7 | 15.8 |
| other socket (cpu 0,64) | 1.20 | 7.2 | 10.7 | 11.7 |
| other node | 2.79 | 15.1 | 2.7 | 11.2 |

Intra-node latency of 1.2-1.4 us, with bandwidth that saturates near the 100 Gb/s link rate, means that intra-node messages use the HCA (rc loopback), not shared memory. This is why the transport became an experiment dimension.

**2026-11-08, jobs 12258964 and 12259215 (pingpong with shared memory; data in `results/pingpong_v1.csv` and `results/diag_ucx.csv`).** `UCX_TLS=sm,self,rc` alone moves only eager messages (up to 8 KB) to shared memory: 0.11 us within a CCD instead of 1.36 us. Above 8 KB, `UCX_PROTO_INFO` shows rendezvous over `rc_mlx5`, so bandwidth stays about 8 GB/s. Adding `UCX_RNDV_THRESH=intra:inf` makes intra-node traffic eager:

| Pair, one-way | 32 KB | 512 KB | 4 MB |
|---|---|---|---|
| same CCD, rc | 4.7 GB/s | 7.7 GB/s | 8.0 GB/s |
| same CCD, shm (intra eager) | 18.4 GB/s | 23.0 GB/s | 26.1 GB/s |
| other socket, rc | 2.2 GB/s | 10.4 GB/s | 11.5 GB/s |
| other socket, shm (intra eager) | 1.1 GB/s | 3.2 GB/s | 3.3 GB/s |

Within a socket, shared memory is about 3x faster than the HCA. Across sockets, it is about 3x slower. So in `shm` mode the cheap level is the socket, not the node. Inter-node bandwidth between 128 KB and 4 MB was about 2.7 GB/s in most runs and about 10 GB/s in one; 8 MB and above reach about 11 GB/s.

**2026-11-08, jobs 12259393 and 12259396 (smoke and pingpong with the final `shm` setting).** Smoke 98/98. Placement check (N = 240, 4x4x3): 0.52 ms per step, `wait` 0.17 ms, compared with 0.87 / 0.64 ms under `rc`. Pingpong, large-message bandwidth and 8 B latency:

| Level | shm | rc | 8 B latency, shm / rc |
|---|---|---|---|
| same CCD | about 25 GB/s | 8-15 GB/s | 0.11 / 1.37 us |
| same socket, other CCD | about 6 GB/s | about 8 GB/s | 0.33 / 1.37 us |
| other socket | about 3 GB/s | about 11 GB/s | 0.78 / 1.18 us |
| other node | about 2.7 GB/s (128 KB-4 MB) | same | 2.75 / 2.75 us |

Under `shm` the cost rises at each level (CCD, socket, node); `rc` makes all intra-node levels look alike. Inter-node bandwidth at 128 KB-4 MB was about 2.7 GB/s in most runs, so use that for the model. rc intra-node numbers vary between runs.

**2026-11-08, job 12259398 (S, scaling, `shm`).** 165 runs, all PASS. Plots: `figs/S_*.png` (`python3 analysis/plot_S.py`).
- Weak scaling, local about 240^3 (does not fit in L3): time per step goes from 9.1 ms at P = 1 to 19.7 ms at P = 8, then stays at about 22 ms up to P = 48. P = 1..8 fills the 8 CCDs of socket 0, so the loss is memory-bandwidth contention within a socket; adding nodes costs almost nothing.
- Weak scaling, local about 60^3 (fits in L3): 0.14 ms to 0.27 ms, a slow rise from communication.
- Strong scaling efficiency at P = 48: 0.49 (N = 240), 0.44 (480), 0.53 (960). The curves do not collapse onto one efficiency vs n/P curve; for N = 480 and 960 most of the drop happens at P = 6..8, inside one socket. n/P where efficiency falls below 0.8: about 4.6e6 (N = 240), 3.2e7 (480), 1.5e8 (960).
- N = 480 at P = 48 (about 37 MB per rank, near the 32 MB L3): efficiency recovers from 0.36 at P = 32 to 0.44.

**2026-11-08, job 12259401 (A, shape, `shm`).** 405 runs, all PASS. Plot: `figs/A_shape.png` (`PYTHONPATH=analysis python3 analysis/plot_A.py`).
- Best grids are near-cubic (2x4x6, 4x3x4, 3x4x4) at every N. Time per step grows with halo volume, and splitting x (the unit-stride direction) adds a penalty on top: at the same halo volume, Px >= 6 grids are the slowest.
- Same surface, different orientation: 48x1x1 is 2.0x (N = 240), 1.8x (480) and 1.3x (960) slower than 1x1x48. Its pack + unpack alone is 32 ms per step at N = 960, against 7.6 ms for 1x1x48.
- Worst vs best: 4.0x (N = 240), 3.4x (480), 2.1x (960). The prediction that shape barely matters at N = 960 was wrong: the slabs are 1.6-2.1x slower even there.
- At N = 480 the fastest grids do not split z (8x6x1, 4x12x1), even though their halo is larger than that of 4x4x3.

**2026-11-08, job 12259403 (B, bad P, `shm`).** 78 runs, all PASS. Plot: `figs/B_badP.png`.
- Every prime P (1x1xP) is slower than every non-prime P from 38 to 48, and slower than or equal to P = 36 (3x3x4).
- N = 240: P = 47 takes 0.71 ms against 0.38 ms at P = 40 (2x4x5) or 48. Using 40 ranks instead of 47 is 1.9x faster.
- N = 960: P = 47 takes 45.5 ms against 30.0 ms at P = 45 (3x3x5). The prediction that the prime penalty almost disappears at N = 960 (about 2%) was wrong; it is about 1.5x, from both `wait` (17 vs 6 ms) and compute (33 vs 26 ms; the slab orientation effect from A).

**2026-11-08, job 12259404 (C, placement x transport).** 108 runs, all PASS. Plot: `figs/C_placement.png`.
- `shm`: placement matters where the extra surface is. With 1x1x48 at N = 240, cyclic is 2.1x slower than block (1.43 vs 0.68 ms); with 4x4x3 it is 1.3x. At N = 960, 4x4x3 does not care (29.2 vs 29.6 ms), and 1x1x48 is 1.2x slower under cyclic.
- `rc`: block and cyclic are close for the slabs (1.34 vs 1.43 ms), so the hierarchy mostly disappears when intra-node traffic goes through the HCA. For 4x4x3 at N = 240, cyclic is even faster than block (0.52 vs 0.85 ms).
- Transport, block placement, N = 240: `shm` is about 2x faster than `rc` (4x4x3 0.40 vs 0.85 ms; 1x1x48 0.68 vs 1.34 ms). At N = 960 the difference is small except for the slabs.
- Packed (8 ranks per CCD): 2.2-4.4x slower at N = 240 and about 3x slower at N = 960, mostly compute (84 vs 25 ms for 4x4x3), because 8 ranks share one L3 and one CCD's memory bandwidth. This is why the baseline is undersubscribed.
- Even with block placement under `shm`, 1x1x48 is 1.7x slower than 4x4x3 at N = 240, although its extra surface stays inside the node: part of it crosses sockets, where shared memory gives only about 3 GB/s.

**2026-11-08, job 12259405 (D, overlap x transport).** 96 runs, all PASS. Plot: `figs/D_overlap.png`.
- `blocking` (one Sendrecv per direction, in sequence) is the slowest at N = 240: 1.7-1.9x slower than `nonblocking` under `shm`.
- N = 960, the clearest case of the progress effect: `overlap` alone barely reduces `wait` (rc 4x4x3: 5.7 -> 7.0 ms; rc 1x1x48: 22.5 -> 19.5 ms), while `overlap_test` cuts it to 0.4 ms and 4.5 ms. Without MPI calls during the interior compute, the rendezvous does not progress.
- The overlap gain is mostly eaten by the separate boundary pass: `comp_bd` costs 4-5 ms at N = 960 (about 18% of the compute for 4x4x3, for about 2% of the cells). Net gain in `t_step`: none for 4x4x3, 15% for rc 1x1x48 (50.4 -> 42.8 ms).
- N = 240: under `rc`, 4x4x3 gains a lot (0.84 -> 0.46 ms with `overlap`, 0.33 ms with `overlap_test`, 2.5x). Under `shm` (intra-node eager) the modes are close. 1x1x48 does not gain in any mode, because its communication (about 1 ms) is much longer than the 0.1-0.2 ms of interior compute it could hide behind.
- Under `shm`, `overlap_test` at N = 960 helps 4x4x3 `wait` (5.1 -> 2.5 ms) but not 1x1x48 (12.9 -> 12.1 ms).

**A model (`analysis/model_A.py`, plot `figs/A_model.png`).** `t_step = a + b*halo + c*log2(Px)`, fitted per N over the 45 grids:

| N | b (ms/MB) | effective GB/s per rank | c (ms per halving of x) | R^2 | R^2, halo only |
|---|---|---|---|---|---|
| 240 | 0.79 | 1.27 | 0.08 | 0.70 | 0.41 |
| 480 | 0.91 | 1.10 | 0.25 | 0.59 | 0.44 |
| 960 | 1.54 | 0.65 | 1.89 | 0.90 | 0.65 |

- Halo volume plus the x-split penalty explains most of the shape effect. The effective bandwidth with all 48 ranks exchanging at once (0.65-1.3 GB/s per rank) is well below the two-point pingpong numbers (3-25 GB/s intra-node, 2.7 GB/s inter-node), so a two-point alpha-beta model underestimates the communication time.
- Off-node bytes per node hardly change across the 45 grids (0.61-0.81 MB at N = 240) while intra-node bytes vary 8x. With block placement, inter-node traffic is set by the node-level decomposition (Section 3), not by the rank grid.

**2026-11-08, job 12259406 (E, packing, `shm`).** 48 runs, all PASS. Plot: `figs/E_packing.png`.
- N = 240: derived datatypes are equal or slower for every face direction (y-faces +21%, z-faces +18%, 4x4x3 +5%, x-faces equal). Even fully contiguous z-faces do not gain at this size.
- N = 960: datatypes are faster for the slabs (y-faces -13%, z-faces -11%, x-faces -3%) and equal for 4x4x3. Manual packing of a contiguous z-face still copies about 7.6 ms per step, which the datatype path avoids.
- Strided x-faces are the expensive case either way: with manual packing, pack + unpack is 32 ms of the 59 ms step at N = 960, more than the communication itself.
- Against the predictions: z-faces gain from datatypes only at large N; x-faces are expensive as expected, but datatypes are not slower than manual packing there.

## 8b. Status and next steps (update as things finish)

- Done: smoke, pingpong, S (12259398), A (12259401), B (12259403), C (12259404), D (12259405), E (12259406). Plots: `figs/S_*.png`, `figs/A_shape.png`, `figs/B_badP.png`, `figs/C_placement.png`, `figs/D_overlap.png`, `figs/E_packing.png`, `figs/A_model.png`.
- Queued: S2 12260302.
- S2 (`scripts/S2_socket.sbatch`): same as S for P <= 16 with the `ALT` map (ranks alternating between sockets). Compare with S to show socket memory-bandwidth contention (tags `strong_alt`, `weak240_alt` in `results/S.csv`).
- After each job: summarize the trend here (8a), make the plot in `analysis/`, commit and push.
- Remaining analysis: S2 vs S comparison plot.
- Decisions made: P = 48 on 3 nodes; N = 240/480/960 in A; `shm` = intra-node eager; pingpong, C and D under both `shm` and `rc`; keep mode name `nonblocking` (the no-overlap control).

## 9. Code and usage

```
src/halo3d.c          the benchmark (one file)
src/pingpong.c        alpha-beta measurement between two ranks
scripts/common.sh     modules, CPU maps, helpers (sourced by every job)
scripts/smoke.sh      correctness sweep (local, or inside an allocation)
scripts/*.sbatch      smoke test, pingpong, experiments S, A-E
scripts/submit_all.sh submit smoke, then everything else once smoke passes
results/              CSV output (one file per experiment) and Slurm logs
analysis/             plotting scripts (run from the repo root)
figs/                 generated plots
```

### Build and run

```bash
# Campus Cluster login node
module load gcc openmpi            # gcc/13.3.0 + openmpi/5.0.1-gcc-13.3.0
make ARCH=-march=znver3
scripts/submit_all.sh              # extra sbatch options are passed to every job

# Local check (any machine with MPI). On the login node there is no InfiniBand,
# so UCX fails to load; use shared memory instead:
OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader scripts/smoke.sh
```

`results/logs/halo-smoke-*.out` should show `smoke: 98/98 passed`, the `lscpu` table (confirm that cores 8c..8c+7 share one L3; otherwise fix `SPREAD`/`PACKED` in `scripts/common.sh`), and a placement listing where each rank sits on a distinct CCD.

### halo3d options

| Option | Meaning | Default |
|---|---|---|
| `--n N` | global cells per dimension (N^3 total) | 256 |
| `--grid PxxPyxPz` | process grid, or `auto` (`MPI_Dims_create`, largest count on z) | auto |
| `--steps S`, `--warmup W` | timed and untimed time steps | 100, 5 |
| `--mode` | `blocking`, `nonblocking`, `overlap`, `overlap_test` | nonblocking |
| `--pack` | `manual`, `datatype` | manual |
| `--test-every T` | `overlap_test`: `MPI_Testall` every T interior z-planes | 1 |
| `--ccd-size C` | cores per CCD, for the placement statistics | 8 |
| `--csv FILE`, `--tag STR` | append one row to FILE | - |
| `--verbose` | print each rank's coordinates, local size, host and core | off |

The program exits with status 1 if the correctness check fails.

### CSV columns

- `tls`: value of `UCX_TLS` for the run, with commas written as `+` (`sm+self+rc` or `rc`).
- `t_step`: max over ranks of the time per step (seconds). This is the main metric.
- `<phase>_avg`, `<phase>_max`: per-step time of each phase, averaged and maximized over ranks. Phases: `pack`, `post` (posting Isend/Irecv), `wait` (`Waitall`, or `Sendrecv` in blocking mode), `unpack`, `comp_int` (interior compute in overlap modes, the whole domain otherwise), `comp_bd` (boundary shell, overlap modes only).
- `halo_bytes_max`, `msgs_max`: largest per-rank halo volume and neighbor count per step.
- `bytes_offnode`, `bytes_samenode`, `bytes_sameccd`: total halo bytes per step summed over all ranks, classified by where the receiving neighbor runs (from `sched_getcpu`, so only meaningful when ranks are pinned).
- `cells_min`, `cells_max`: smallest and largest subdomain (load imbalance when N is not divisible).
- `rel_err`: max relative error against the exact discrete solution (PASS below 1e-10).

`results/pingpong.csv`: `tag` (placement), `tls`, `bytes`, `reps` round trips per trial, `trials`, and the one-way time `t_min` / `t_med` (min and median over trials, seconds) with the matching `gbps_min` / `gbps_med`. Fit alpha and beta to `t_min` or `t_med`.
