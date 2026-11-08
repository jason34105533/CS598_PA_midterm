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

Inter-node traffic is identical, because the node-level decomposition is 1x1x3 in both cases. The extra surface of the slab is all intra-node (shared memory). This is why 4x4x3 is the "good" grid in C, D and E rather than 3x4x4 (same surface, but a z-layer of 12 ranks does not line up with a 16-rank node). With cyclic rank distribution, the same 1x1x48 sends every face off-node, which is the worst case. So decomposition quality has to be judged at each level: node, CCD, and core. The CCD level comes in naturally here.

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
| C | Same grid, different placement: block vs cyclic across nodes; inside a node, one rank per CCD vs packed into few CCDs. | Does it matter where the extra surface lands (node / CCD / core)? |
| D | blocking vs overlap vs overlap_test, on the best and worst decomposition from A. | How much communication can overlap hide? Does it really overlap (compare Waitall time with and without `MPI_Test`)? |
| E | manual packing vs derived datatypes, broken down by face direction. | How expensive are non-contiguous faces? Does the answer depend on the decomposition? |
| F (stretch) | Irregular decomposition for prime P. | Can it beat 1x1xP and beat dropping ranks? |

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

## 9. Code and usage

```
src/halo3d.c          the benchmark (one file)
src/pingpong.c        alpha-beta measurement between two ranks
scripts/common.sh     modules, CPU maps, helpers (sourced by every job)
scripts/smoke.sh      correctness sweep (local, or inside an allocation)
scripts/*.sbatch      smoke test, pingpong, experiments S, A-E
scripts/submit_all.sh submit smoke, then everything else once smoke passes
results/              CSV output (one file per experiment) and Slurm logs
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

- `t_step`: max over ranks of the time per step (seconds). This is the main metric.
- `<phase>_avg`, `<phase>_max`: per-step time of each phase, averaged and maximized over ranks. Phases: `pack`, `post` (posting Isend/Irecv), `wait` (`Waitall`, or `Sendrecv` in blocking mode), `unpack`, `comp_int` (interior compute in overlap modes, the whole domain otherwise), `comp_bd` (boundary shell, overlap modes only).
- `halo_bytes_max`, `msgs_max`: largest per-rank halo volume and neighbor count per step.
- `bytes_offnode`, `bytes_samenode`, `bytes_sameccd`: total halo bytes per step summed over all ranks, classified by where the receiving neighbor runs (from `sched_getcpu`, so only meaningful when ranks are pinned).
- `cells_min`, `cells_max`: smallest and largest subdomain (load imbalance when N is not divisible).
- `rel_err`: max relative error against the exact discrete solution (PASS below 1e-10).
