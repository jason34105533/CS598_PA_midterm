/*
 * halo3d_irr.c -- the halo3d benchmark (same stencil, same correctness check) on a
 * layered, irregular decomposition, for P with no good 3D factorization.
 *
 * --layers PxxPy,PxxPy,...  splits z into layers, one per entry; layer l is split
 * into a Px x Py grid in xy. Layer thickness is proportional to its rank count, so
 * every rank gets about N^3/P cells. Example: P = 47 = 16 + 16 + 15 with
 * --layers 4x4,4x4,3x5. Across a layer boundary the xy grids do not line up, so a
 * rank can have several z-neighbors; neighbors and exchanged regions are found by
 * intersecting boxes. A regular grid is the special case of equal layers.
 *
 * Only the nonblocking exchange with manual packing is implemented (as in B).
 * Rank numbering: layer by layer, x fastest within a layer.
 */
#define _GNU_SOURCE
#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef __linux__
#include <sched.h>
#endif

#define MAXL 64
#define MAXNB 64

enum { T_PACK, T_POST, T_WAIT, T_UNPACK, T_COMP, T_NUM };
static const char *tname[T_NUM] = {"pack", "post", "wait", "unpack", "comp_int"};

typedef struct { int off[3], n[3]; } Box;

typedef struct {
    int peer;           /* neighbor rank */
    int lo[3], hi[3];   /* local index range sent (1-based, inclusive) */
    int rlo[3], rhi[3]; /* local index range received into (ghost cells) */
    int count;
    double *sbuf, *rbuf;
} Msg;

static int rank_, size_;
static int N = 240, steps = 100, warmup = 5, ccd_size = 8, verbose = 0;
static double r = 0.125;
static const char *csv = NULL, *tag = "-", *layers_arg = NULL;
static int L, lpx[MAXL], lpy[MAXL];

static void die(const char *msg)
{
    if (rank_ == 0) fprintf(stderr, "error: %s\n", msg);
    MPI_Abort(MPI_COMM_WORLD, 1);
}

static void split(int n, int p, int c, int *cnt, int *off)
{
    int base = n / p, rem = n % p;
    *cnt = base + (c < rem);
    *off = c * base + (c < rem ? c : rem);
}

/* Box of rank q: find its layer, then its xy position in the layer grid. */
static Box box_of(int q)
{
    int base = 0, l = 0;
    while (q >= base + lpx[l] * lpy[l]) base += lpx[l] * lpy[l], l++;
    int before = base, upto = base + lpx[l] * lpy[l];
    Box b;
    /* z: layer thickness proportional to its rank count */
    int z0 = (int)((long)N * before / size_), z1 = (int)((long)N * upto / size_);
    b.off[2] = z0; b.n[2] = z1 - z0;
    int i = q - base, cx = i % lpx[l], cy = i / lpx[l];
    split(N, lpx[l], cx, &b.n[0], &b.off[0]);
    split(N, lpy[l], cy, &b.n[1], &b.off[1]);
    return b;
}

#define IDX(i, j, k) ((size_t)(i) + sx * (size_t)(j) + sy * (size_t)(k))
static size_t sx, sy, len;

/* Copy the local index box [lo, hi] between u and a contiguous buffer. */
static void region_copy(double *u, const int *lo, const int *hi, double *buf, int to_buf)
{
    size_t m = 0;
    int w = hi[0] - lo[0] + 1;
    for (int k = lo[2]; k <= hi[2]; k++)
        for (int j = lo[1]; j <= hi[1]; j++) {
            double *row = u + IDX(lo[0], j, k);
            if (w == 1) { if (to_buf) buf[m] = *row; else *row = buf[m]; }
            else if (to_buf) memcpy(buf + m, row, sizeof(double) * w);
            else memcpy(row, buf + m, sizeof(double) * w);
            m += w;
        }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &size_);

    for (int a = 1; a < argc; a++) {
        const char *k = argv[a], *v = a + 1 < argc ? argv[a + 1] : NULL;
        if (!strcmp(k, "--verbose")) { verbose = 1; continue; }
        if (!v) die("missing option value");
        a++;
        if (!strcmp(k, "--n")) N = atoi(v);
        else if (!strcmp(k, "--layers")) layers_arg = v;
        else if (!strcmp(k, "--steps")) steps = atoi(v);
        else if (!strcmp(k, "--warmup")) warmup = atoi(v);
        else if (!strcmp(k, "--ccd-size")) ccd_size = atoi(v);
        else if (!strcmp(k, "--csv")) csv = v;
        else if (!strcmp(k, "--tag")) tag = v;
        else die("unknown option (--n --layers --steps --warmup --ccd-size --csv --tag --verbose)");
    }
    if (!layers_arg) die("--layers PxxPy,PxxPy,... is required");
    {
        char tmp[512];
        snprintf(tmp, sizeof tmp, "%s", layers_arg);
        int total = 0;
        for (char *s = strtok(tmp, ","); s; s = strtok(NULL, ",")) {
            if (L == MAXL || sscanf(s, "%dx%d", &lpx[L], &lpy[L]) != 2 || lpx[L] < 1 || lpy[L] < 1)
                die("bad --layers entry");
            total += lpx[L] * lpy[L];
            L++;
        }
        if (total != size_) die("sum of Px*Py over layers must equal the number of ranks");
    }

    Box me = box_of(rank_);
    for (int d = 0; d < 3; d++) if (me.n[d] < 1) die("a subdomain is empty; N too small");
    sx = (size_t)me.n[0] + 2;
    sy = sx * ((size_t)me.n[1] + 2);
    len = sy * ((size_t)me.n[2] + 2);

    /* Neighbors: every box that shares part of a face with mine. */
    Msg msg[MAXNB];
    int nm = 0;
    for (int q = 0; q < size_; q++) {
        if (q == rank_) continue;
        Box o = box_of(q);
        for (int d = 0; d < 3; d++) {
            int plus;
            if (o.off[d] == me.off[d] + me.n[d]) plus = 1;
            else if (o.off[d] + o.n[d] == me.off[d]) plus = 0;
            else continue;
            int lo[3], hi[3], ok = 1;
            for (int e = 0; e < 3; e++) {
                if (e == d) continue;
                int a0 = me.off[e] > o.off[e] ? me.off[e] : o.off[e];
                int a1 = me.off[e] + me.n[e] < o.off[e] + o.n[e] ? me.off[e] + me.n[e] : o.off[e] + o.n[e];
                if (a1 <= a0) { ok = 0; break; }
                lo[e] = a0 - me.off[e] + 1;
                hi[e] = a1 - me.off[e];
            }
            if (!ok) continue;
            if (nm == MAXNB) die("too many neighbors");
            Msg *m = &msg[nm++];
            m->peer = q;
            memcpy(m->lo, lo, sizeof lo); memcpy(m->hi, hi, sizeof hi);
            memcpy(m->rlo, lo, sizeof lo); memcpy(m->rhi, hi, sizeof hi);
            m->lo[d] = m->hi[d] = plus ? me.n[d] : 1;
            m->rlo[d] = m->rhi[d] = plus ? me.n[d] + 1 : 0;
            m->count = 1;
            for (int e = 0; e < 3; e++) m->count *= m->hi[e] - m->lo[e] + 1;
            m->sbuf = malloc(sizeof(double) * m->count);
            m->rbuf = malloc(sizeof(double) * m->count);
            if (!m->sbuf || !m->rbuf) die("out of memory (buffers)");
        }
    }

    double *u = calloc(len, sizeof(double)), *un = calloc(len, sizeof(double));
    if (!u || !un) die("out of memory (field)");
    const double h = 1.0 / (N + 1);
    for (int k = 1; k <= me.n[2]; k++)
        for (int j = 1; j <= me.n[1]; j++)
            for (int i = 1; i <= me.n[0]; i++)
                u[IDX(i, j, k)] = sin(M_PI * h * (me.off[0] + i)) * sin(M_PI * h * (me.off[1] + j)) *
                                  sin(M_PI * h * (me.off[2] + k));

    /* Placement statistics, as in halo3d. */
    int *node_of = malloc(sizeof(int) * size_), *cpu_of = malloc(sizeof(int) * size_), nodes;
    {
        MPI_Comm nc;
        int lr, leader = rank_, cpu = -1;
        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &nc);
        MPI_Comm_rank(nc, &lr);
        MPI_Allreduce(MPI_IN_PLACE, &leader, 1, MPI_INT, MPI_MIN, nc);
        MPI_Comm_free(&nc);
#ifdef __linux__
        cpu = sched_getcpu();
#endif
        MPI_Allgather(&leader, 1, MPI_INT, node_of, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Allgather(&cpu, 1, MPI_INT, cpu_of, 1, MPI_INT, MPI_COMM_WORLD);
        int is_leader = lr == 0;
        MPI_Allreduce(&is_leader, &nodes, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }
    double bytes[3] = {0, 0, 0};
    for (int i = 0; i < nm; i++) {
        int q = msg[i].peer;
        double b = 8.0 * msg[i].count;
        if (node_of[q] != node_of[rank_]) bytes[0] += b;
        else if (ccd_size > 0 && cpu_of[q] >= 0 && cpu_of[q] / ccd_size == cpu_of[rank_] / ccd_size) bytes[2] += b;
        else bytes[1] += b;
    }
    if (verbose) {
        char line[200];
        snprintf(line, sizeof line, "rank %4d  box off (%d,%d,%d) n %dx%dx%d  neighbors %d  cpu %d node %d",
                 rank_, me.off[0], me.off[1], me.off[2], me.n[0], me.n[1], me.n[2], nm, cpu_of[rank_], node_of[rank_]);
        char *all = rank_ == 0 ? malloc(sizeof line * size_) : NULL;
        MPI_Gather(line, sizeof line, MPI_CHAR, all, sizeof line, MPI_CHAR, 0, MPI_COMM_WORLD);
        if (rank_ == 0) { for (int p = 0; p < size_; p++) printf("%s\n", all + p * sizeof line); free(all); }
    }

    MPI_Request req[2 * MAXNB];
    double t[T_NUM] = {0}, ttot = 0;
    for (int s = -warmup; s < steps; s++) {
        if (s == 0) { memset(t, 0, sizeof t); MPI_Barrier(MPI_COMM_WORLD); ttot = MPI_Wtime(); }
        double t0 = MPI_Wtime(), t1;
#define LAP(p) (t1 = MPI_Wtime(), t[p] += t1 - t0, t0 = t1)
        for (int i = 0; i < nm; i++) region_copy(u, msg[i].lo, msg[i].hi, msg[i].sbuf, 1);
        LAP(T_PACK);
        /* one message per neighbor pair and direction, so tag 0 is enough */
        for (int i = 0; i < nm; i++)
            MPI_Irecv(msg[i].rbuf, msg[i].count, MPI_DOUBLE, msg[i].peer, 0, MPI_COMM_WORLD, &req[i]);
        for (int i = 0; i < nm; i++)
            MPI_Isend(msg[i].sbuf, msg[i].count, MPI_DOUBLE, msg[i].peer, 0, MPI_COMM_WORLD, &req[nm + i]);
        LAP(T_POST);
        MPI_Waitall(2 * nm, req, MPI_STATUSES_IGNORE);
        LAP(T_WAIT);
        for (int i = 0; i < nm; i++) region_copy(u, msg[i].rlo, msg[i].rhi, msg[i].rbuf, 0);
        LAP(T_UNPACK);
        for (int k = 1; k <= me.n[2]; k++)
            for (int j = 1; j <= me.n[1]; j++) {
                const size_t base = sx * (size_t)j + sy * (size_t)k;
                for (int i = 1; i <= me.n[0]; i++) {
                    const size_t p = base + (size_t)i;
                    un[p] = u[p] + r * (u[p - 1] + u[p + 1] + u[p - sx] + u[p + sx] + u[p - sy] +
                                        u[p + sy] - 6.0 * u[p]);
                }
            }
        LAP(T_COMP);
#undef LAP
        double *tmp = u; u = un; un = tmp;
    }
    ttot = (MPI_Wtime() - ttot) / steps;

    /* Correctness: exact discrete solution g^K u0. */
    const double s2 = sin(M_PI * h / 2), G = pow(1.0 - 12.0 * r * s2 * s2, warmup + steps);
    double err = 0, gerr;
    for (int k = 1; k <= me.n[2]; k++)
        for (int j = 1; j <= me.n[1]; j++)
            for (int i = 1; i <= me.n[0]; i++) {
                double ex = G * sin(M_PI * h * (me.off[0] + i)) * sin(M_PI * h * (me.off[1] + j)) *
                            sin(M_PI * h * (me.off[2] + k));
                double e = fabs(u[IDX(i, j, k)] - ex);
                if (e > err) err = e;
            }
    MPI_Allreduce(&err, &gerr, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    double rel_err = gerr / G;
    int pass = rel_err < 1e-10;

    double tavg[T_NUM], tmax[T_NUM], tstep;
    for (int p = 0; p < T_NUM; p++) t[p] /= steps;
    MPI_Reduce(t, tavg, T_NUM, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(t, tmax, T_NUM, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ttot, &tstep, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    double my_bytes = bytes[0] + bytes[1] + bytes[2], max_bytes, sum_bytes[3];
    long my_cells = (long)me.n[0] * me.n[1] * me.n[2], max_cells, min_cells;
    int max_msgs;
    MPI_Reduce(&my_bytes, &max_bytes, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(bytes, sum_bytes, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&nm, &max_msgs, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&my_cells, &max_cells, 1, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&my_cells, &min_cells, 1, MPI_LONG, MPI_MIN, 0, MPI_COMM_WORLD);

    if (rank_ == 0) {
        for (int p = 0; p < T_NUM; p++) tavg[p] /= size_;
        printf("halo3d_irr  N=%d  layers=%s  P=%d  nodes=%d  steps=%d\n", N, layers_arg, size_, nodes, steps);
        printf("  time/step %.3f ms   (avg | max over ranks, ms)\n", 1e3 * tstep);
        for (int p = 0; p < T_NUM; p++) printf("    %-9s %9.4f | %9.4f\n", tname[p], 1e3 * tavg[p], 1e3 * tmax[p]);
        printf("  halo bytes/step: max per rank %.0f (max %d msgs); total off-node %.0f, same-node %.0f, same-CCD %.0f\n",
               max_bytes, max_msgs, sum_bytes[0], sum_bytes[1], sum_bytes[2]);
        printf("  local cells min %ld max %ld   rel_err %.3e  %s\n", min_cells, max_cells, rel_err, pass ? "PASS" : "FAIL");
        if (csv) {
            FILE *fp = fopen(csv, "a");
            if (fp) {
                fseek(fp, 0, SEEK_END);
                if (ftell(fp) == 0) {
                    fprintf(fp, "tag,tls,jobid,layers,N,P,nodes,steps,t_step");
                    for (int p = 0; p < T_NUM; p++) fprintf(fp, ",%s_avg,%s_max", tname[p], tname[p]);
                    fprintf(fp, ",halo_bytes_max,msgs_max,bytes_offnode,bytes_samenode,bytes_sameccd,"
                                "cells_min,cells_max,rel_err\n");
                }
                const char *job = getenv("SLURM_JOB_ID");
                char tls[128], lay[512];
                snprintf(tls, sizeof tls, "%s", getenv("UCX_TLS") ? getenv("UCX_TLS") : "-");
                for (char *c = tls; *c; c++) if (*c == ',') *c = '+';
                snprintf(lay, sizeof lay, "%s", layers_arg);
                for (char *c = lay; *c; c++) if (*c == ',') *c = '+';
                fprintf(fp, "%s,%s,%s,%s,%d,%d,%d,%d,%.6e", tag, tls, job ? job : "-", lay, N, size_, nodes, steps, tstep);
                for (int p = 0; p < T_NUM; p++) fprintf(fp, ",%.6e,%.6e", tavg[p], tmax[p]);
                fprintf(fp, ",%.0f,%d,%.0f,%.0f,%.0f,%ld,%ld,%.3e\n", max_bytes, max_msgs, sum_bytes[0],
                        sum_bytes[1], sum_bytes[2], min_cells, max_cells, rel_err);
                fclose(fp);
            }
        }
    }
    for (int i = 0; i < nm; i++) { free(msg[i].sbuf); free(msg[i].rbuf); }
    free(u); free(un); free(node_of); free(cpu_of);
    MPI_Finalize();
    return pass ? 0 : 1;
}
