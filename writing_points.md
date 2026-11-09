# Writing points

Explanations that worked well, kept for the report and slides. Numbers come from README 8a.

## Storyline

Suggested order: S/S2 → A (with the model, E as its mechanism) → B/F → D. C and E get one figure each.

## Transport: the hierarchy only exists with shared memory

- The site default `UCX_TLS=rc` sends intra-node messages through the HCA, so same CCD, other CCD and other socket all look alike (about 1.2-1.4 us, 8-11 GB/s).
- With `sm` and intra-node eager, the hierarchy appears: same CCD about 25 GB/s, same socket about 6 GB/s, other socket about 3 GB/s, other node about 2.7 GB/s.
- This setting is not the default, which is a finding in itself (worth one slide).

## S/S2: the socket's memory bandwidth limits large problems

- Weak scaling at 240^3 per rank doubles the time per step from P = 1 to 8 (filling the 8 CCDs of socket 0), then stays flat across nodes. Adding nodes costs almost nothing.
- S2 places the same ranks alternately on the two sockets: efficiency at P = 6-12 rises by up to 0.3, and the curves meet at P = 16. This confirms the S result directly.

## A + E: shape, and why splitting x costs extra

- Near-cubic grids are best. Worst vs best is 2-4x.
- Model `t = a + b*halo + c*log2(Px)` explains most of it (R^2 0.59-0.90). The effective bandwidth with 48 ranks exchanging at once is 0.65-1.3 GB/s per rank, well below the two-point pingpong numbers.
- E gives the mechanism for the `c*log2(Px)` term: the x-faces are strided, and 48x1x1 spends 32 ms per step in pack + unpack against 7.6 ms for 1x1x48 (N = 960). Datatypes vs manual packing differ by only 10-20%; present E as support for A, not as a topic of its own.

## B/F: the predictions that were wrong are the best points

- We predicted the prime-P penalty would almost vanish at N = 960. It is 1.5x (P = 47 vs 45).
- F: a layered irregular decomposition for prime P beats 1x1xP by about 1.5x at N = 960 and ties with dropping to fewer ranks. When communication dominates (N = 240), dropping to a P with a good factorization is simpler and faster.

## D: overlap needs three conditions

Overlap can save at most min(communication, compute), so it only pays when the two are of similar size.

| Case | compute | wait | what overlap could save |
|---|---|---|---|
| N = 960, 4x4x3 | about 24 ms | about 5 ms | at most 17% of the step (compute is too **large**) |
| N = 240, 1x1x48, rc | 0.18 ms | 1 ms | at most 0.18 ms (here compute is too small) |
| N = 240, 4x4x3, rc | 0.18 ms | 0.6 ms | closest to balanced; **actually 2.5x faster** (0.84 → 0.33 ms) |
| N = 480, 4x4x3, rc | 1.65 ms | 2.6 ms | **2.1x faster** (4.22 → 2.04 ms); `wait` 2.6 → 0.2 ms |

Why the gain is small elsewhere:

1. **The ratio is lopsided.** Our N values sit on either side of the balanced region.
2. **Under `shm` there is no background transfer.** Overlap needs hardware that moves the data while the CPU computes. Under rc the HCA does it by DMA. Under shm the transfer is a CPU memcpy: `post` (inside `MPI_Isend`) rises from 0.002 ms under rc to 0.05-0.2 ms at N = 240 and 0.4-1.2 ms at N = 960, and the rest needs the receiver's progress. This is why overlap works only under rc.
3. **The separate boundary pass is expensive.** About 2% of the cells cost about 4 ms at N = 960 (strided x-faces, data read again), which eats the hidden time.
4. **Part of `wait` is waiting for slower neighbors.** shm N = 960 4x4x3: 3.3 ms on average, 5.1 ms on the slowest rank. Overlap cannot hide that.

Gain vs N under rc (`figs/D2_gain.png`): 4x4x3 2.5x / 2.1x / 1.0x at N = 240 / 480 / 960; 1x1x48 1.03x / 1.32x / 1.18x. Under shm no gain at any N.

The peel order (boundary first, send, then interior, the idea from class) gives the same times as `overlap_test`: both hide the same exchange behind the same interior compute.

Conclusion: overlap pays off when (1) hardware can move data in the background (HCA DMA under rc), (2) MPI is driven during compute (`MPI_Testall`), and (3) communication and compute are of similar size. The one case in our data that meets all three is 2.5x faster.

## Concurrency explains the A model's bandwidth

- Single-pair pingpong between nodes gives 2.4-2.7 GB/s. With 16 ranks per node exchanging off-node at once, each gets 0.76 GB/s, and together 12 GB/s, the 100 Gb/s link.
- The A model fitted 0.65-1.3 GB/s per rank. So its effective bandwidth is the node's link shared by 16 ranks: a two-point alpha-beta model must divide beta's bandwidth by the ranks sharing the link.

## Index lists (indirect addressing)

- Packing through precomputed index lists costs the same as structured loops for large faces, and as an MPI indexed type it performs like a subarray type. It does not make strided x-faces cheaper (the same cache lines are touched), but it gives one code path for every face and for irregular decompositions at no cost at large N.

## C: why the baseline is undersubscribed

- Packed placement (8 ranks per CCD) is 2-4x slower, mostly in compute, because 8 ranks share one L3 and one CCD's memory bandwidth. This justifies using 16 of 128 cores per node.

## State once, without elaborating

- 16 of 128 cores per node, one rank per CCD (justified by C).
- Each configuration is run 3 times; we report the median.
- `shm` mode is not the site default (`rc`).
