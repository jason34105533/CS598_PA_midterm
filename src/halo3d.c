/*
 * halo3d.c -- 3D heat diffusion (explicit, 7-point stencil) on the unit cube
 * with a configurable halo exchange, for studying domain decomposition,
 * overlap, data packing and rank placement.
 *
 *   u_new = u + r * (sum of 6 face neighbors - 6 u),  Dirichlet zero boundary
 *
 * Correctness: with u0 = sin(pi x) sin(pi y) sin(pi z) the discrete solution
 * after K steps is exactly g^K u0, g = 1 - 12 r sin^2(pi h / 2).
 *
 * Layout: x is unit stride. Local array has one ghost layer:
 *   idx(i,j,k) = i + sx*j + sy*k,  sx = nx+2,  sy = (nx+2)(ny+2)
 * Rank numbering: rank = cx + Px*(cy + Py*cz)  (x fastest).
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

enum { MODE_BLOCKING, MODE_NONBLOCKING, MODE_OVERLAP, MODE_OVERLAP_TEST };
static const char *mode_name[] = {"blocking", "nonblocking", "overlap", "overlap_test"};

enum { PACK_MANUAL, PACK_DATATYPE };
static const char *pack_name[] = {"manual", "datatype"};

/* comp_int = interior compute in overlap modes, whole domain otherwise */
enum { T_PACK, T_POST, T_WAIT, T_UNPACK, T_CINT, T_CBD, T_NUM };
static const char *tname[T_NUM] = {"pack", "post", "wait", "unpack", "comp_int", "comp_bd"};

typedef struct {
    int N, steps, warmup, mode, pack, test_every, ccd_size, verbose;
    int P[3]; /* 0 = auto */
    double r;
    const char *csv, *tag;
} Opts;

typedef struct {
    int n[3], off[3], c[3];
    int nbr[6];         /* faces: -x +x -y +y -z +z */
    int fcount[6];      /* doubles per face */
    size_t sx, sy, len; /* strides and total array length */
    double *sbuf[6], *rbuf[6];
    MPI_Datatype stype[6], rtype[6];
} Dom;

#define IDX(d, i, j, k) ((size_t)(i) + (d)->sx * (size_t)(j) + (d)->sy * (size_t)(k))

static int rank_, size_;

static void die(const char *msg)
{
    if (rank_ == 0) fprintf(stderr, "error: %s\n", msg);
    MPI_Abort(MPI_COMM_WORLD, 1);
}

static void usage(void)
{
    if (rank_ == 0)
        fprintf(stderr,
                "usage: halo3d [options]\n"
                "  --n N            global cells per dimension (default 256)\n"
                "  --grid PxxPyxPz  process grid, e.g. 4x4x4, or 'auto' (default auto)\n"
                "  --steps S        timed steps (default 100)\n"
                "  --warmup W       untimed warmup steps (default 5)\n"
                "  --mode M         blocking | nonblocking | overlap | overlap_test\n"
                "  --pack K         manual | datatype\n"
                "  --test-every T   overlap_test: MPI_Testall every T interior z-planes (default 1)\n"
                "  --ccd-size C     cores per CCD for placement stats (default 8, 0 = off)\n"
                "  --r R            diffusion number (default 0.125, must be <= 1/6)\n"
                "  --csv FILE       append one result row to FILE\n"
                "  --tag STR        label stored in the csv row\n"
                "  --verbose        print rank placement\n");
}

static int lookup(const char *s, const char **names, int n)
{
    for (int i = 0; i < n; i++)
        if (strcmp(s, names[i]) == 0) return i;
    return -1;
}

static void parse(int argc, char **argv, Opts *o)
{
    *o = (Opts){.N = 256, .steps = 100, .warmup = 5, .mode = MODE_NONBLOCKING,
                .pack = PACK_MANUAL, .test_every = 1, .ccd_size = 8, .r = 0.125,
                .csv = NULL, .tag = "-"};
    for (int a = 1; a < argc; a++) {
        const char *k = argv[a];
        const char *v = (a + 1 < argc) ? argv[a + 1] : NULL;
        if (!strcmp(k, "--verbose")) { o->verbose = 1; continue; }
        if (!strcmp(k, "--help")) { usage(); MPI_Finalize(); exit(0); }
        if (!v) { usage(); die("missing option value"); }
        a++;
        if (!strcmp(k, "--n")) o->N = atoi(v);
        else if (!strcmp(k, "--steps")) o->steps = atoi(v);
        else if (!strcmp(k, "--warmup")) o->warmup = atoi(v);
        else if (!strcmp(k, "--test-every")) o->test_every = atoi(v);
        else if (!strcmp(k, "--ccd-size")) o->ccd_size = atoi(v);
        else if (!strcmp(k, "--r")) o->r = atof(v);
        else if (!strcmp(k, "--csv")) o->csv = v;
        else if (!strcmp(k, "--tag")) o->tag = v;
        else if (!strcmp(k, "--mode")) {
            if ((o->mode = lookup(v, mode_name, 4)) < 0) die("unknown --mode");
        } else if (!strcmp(k, "--pack")) {
            if ((o->pack = lookup(v, pack_name, 2)) < 0) die("unknown --pack");
        } else if (!strcmp(k, "--grid")) {
            if (strcmp(v, "auto") && sscanf(v, "%dx%dx%d", &o->P[0], &o->P[1], &o->P[2]) != 3)
                die("bad --grid, expected PxxPyxPz or auto");
        } else { usage(); die("unknown option"); }
    }
    if (o->N < 1 || o->steps < 1 || o->warmup < 0 || o->test_every < 1) die("bad numeric option");
    if (o->r <= 0 || o->r > 1.0 / 6.0) die("--r must be in (0, 1/6]");
}

/* Split N cells over P parts as evenly as possible; part c gets n cells at off. */
static void split(int N, int P, int c, int *n, int *off)
{
    int base = N / P, rem = N % P;
    *n = base + (c < rem);
    *off = c * base + (c < rem ? c : rem);
}

static void setup_domain(Dom *d, const Opts *o)
{
    const int *P = o->P;
    d->c[0] = rank_ % P[0];
    d->c[1] = (rank_ / P[0]) % P[1];
    d->c[2] = rank_ / (P[0] * P[1]);
    for (int dim = 0; dim < 3; dim++) split(o->N, P[dim], d->c[dim], &d->n[dim], &d->off[dim]);

    for (int f = 0; f < 6; f++) {
        int dim = f / 2, step = (f % 2) ? 1 : -1, cc[3] = {d->c[0], d->c[1], d->c[2]};
        cc[dim] += step;
        d->nbr[f] = (cc[dim] < 0 || cc[dim] >= P[dim])
                        ? MPI_PROC_NULL
                        : cc[0] + P[0] * (cc[1] + P[1] * cc[2]);
        d->fcount[f] = d->n[(dim + 1) % 3] * d->n[(dim + 2) % 3];
    }

    d->sx = (size_t)d->n[0] + 2;
    d->sy = d->sx * ((size_t)d->n[1] + 2);
    d->len = d->sy * ((size_t)d->n[2] + 2);

    for (int f = 0; f < 6; f++) {
        d->sbuf[f] = malloc(sizeof(double) * d->fcount[f]);
        d->rbuf[f] = malloc(sizeof(double) * d->fcount[f]);
        if (!d->sbuf[f] || !d->rbuf[f]) die("out of memory (face buffers)");
    }

    /* Subarray types address the face directly inside the local array.
       MPI_ORDER_C wants dimensions slowest first: (z, y, x). */
    int sizes[3] = {d->n[2] + 2, d->n[1] + 2, d->n[0] + 2};
    for (int f = 0; f < 6; f++) {
        int dim = f / 2, cd = 2 - dim, plus = f % 2;
        int sub[3] = {d->n[2], d->n[1], d->n[0]}, st[3] = {1, 1, 1};
        sub[cd] = 1;
        st[cd] = plus ? d->n[dim] : 1; /* layer sent */
        MPI_Type_create_subarray(3, sizes, sub, st, MPI_ORDER_C, MPI_DOUBLE, &d->stype[f]);
        MPI_Type_commit(&d->stype[f]);
        st[cd] = plus ? d->n[dim] + 1 : 0; /* ghost layer received */
        MPI_Type_create_subarray(3, sizes, sub, st, MPI_ORDER_C, MPI_DOUBLE, &d->rtype[f]);
        MPI_Type_commit(&d->rtype[f]);
    }
}

/* Copy face f between the array and a contiguous buffer.
   layer = index along the face normal; to_buf = 1 packs, 0 unpacks. */
static void face_copy(const Dom *d, double *u, int f, int layer, double *buf, int to_buf)
{
    const int nx = d->n[0], ny = d->n[1], nz = d->n[2];
    size_t m = 0;
    switch (f / 2) {
    case 0: /* x-face: strided gather, one element per (j,k) */
        for (int k = 1; k <= nz; k++)
            for (int j = 1; j <= ny; j++, m++) {
                size_t p = IDX(d, layer, j, k);
                if (to_buf) buf[m] = u[p]; else u[p] = buf[m];
            }
        break;
    case 1: /* y-face: one x-row per k */
        for (int k = 1; k <= nz; k++, m += nx) {
            double *row = u + IDX(d, 1, layer, k);
            if (to_buf) memcpy(buf + m, row, sizeof(double) * nx);
            else memcpy(row, buf + m, sizeof(double) * nx);
        }
        break;
    case 2: /* z-face: one x-row per j */
        for (int j = 1; j <= ny; j++, m += nx) {
            double *row = u + IDX(d, 1, j, layer);
            if (to_buf) memcpy(buf + m, row, sizeof(double) * nx);
            else memcpy(row, buf + m, sizeof(double) * nx);
        }
        break;
    }
}

static void pack_all(const Dom *d, double *u)
{
    for (int f = 0; f < 6; f++)
        if (d->nbr[f] != MPI_PROC_NULL)
            face_copy(d, u, f, (f % 2) ? d->n[f / 2] : 1, d->sbuf[f], 1);
}

static void unpack_all(const Dom *d, double *u)
{
    for (int f = 0; f < 6; f++)
        if (d->nbr[f] != MPI_PROC_NULL)
            face_copy(d, u, f, (f % 2) ? d->n[f / 2] + 1 : 0, d->rbuf[f], 0);
}

/* Update all cells with i0<=i<=i1, j0<=j<=j1, k0<=k<=k1 (empty if a range is empty). */
static void stencil_box(const Dom *d, const double *restrict u, double *restrict un, double r,
                        int i0, int i1, int j0, int j1, int k0, int k1)
{
    const size_t sx = d->sx, sy = d->sy;
    for (int k = k0; k <= k1; k++)
        for (int j = j0; j <= j1; j++) {
            const size_t base = sx * (size_t)j + sy * (size_t)k;
            for (int i = i0; i <= i1; i++) {
                const size_t p = base + (size_t)i;
                un[p] = u[p] + r * (u[p - 1] + u[p + 1] + u[p - sx] + u[p + sx] +
                                    u[p - sy] + u[p + sy] - 6.0 * u[p]);
            }
        }
}

/* Cells that do not touch any ghost cell. If reqs != NULL, call MPI_Testall
   every test_every z-planes until the exchange has completed. */
static void compute_interior(const Dom *d, const double *u, double *un, double r,
                             MPI_Request *reqs, int test_every)
{
    const int nx = d->n[0], ny = d->n[1], nz = d->n[2];
    if (!reqs) {
        stencil_box(d, u, un, r, 2, nx - 1, 2, ny - 1, 2, nz - 1);
        return;
    }
    int done = 0;
    for (int k = 2; k <= nz - 1; k++) {
        stencil_box(d, u, un, r, 2, nx - 1, 2, ny - 1, k, k);
        if (!done && (k - 1) % test_every == 0) MPI_Testall(12, reqs, &done, MPI_STATUSES_IGNORE);
    }
}

/* The one-cell shell next to the ghost layer, split into disjoint pieces. */
static void compute_boundary(const Dom *d, const double *u, double *un, double r)
{
    const int nx = d->n[0], ny = d->n[1], nz = d->n[2];
    stencil_box(d, u, un, r, 1, nx, 1, ny, 1, 1);
    if (nz > 1) stencil_box(d, u, un, r, 1, nx, 1, ny, nz, nz);
    stencil_box(d, u, un, r, 1, nx, 1, 1, 2, nz - 1);
    if (ny > 1) stencil_box(d, u, un, r, 1, nx, ny, ny, 2, nz - 1);
    stencil_box(d, u, un, r, 1, 1, 2, ny - 1, 2, nz - 1);
    if (nx > 1) stencil_box(d, u, un, r, nx, nx, 2, ny - 1, 2, nz - 1);
}

/* One time step: exchange halos of u, write u_new into un. Adds phase times to t. */
static void step(const Dom *d, const Opts *o, double *u, double *un, double *t)
{
    const int manual = (o->pack == PACK_MANUAL);
    MPI_Request req[12];
    double t0 = MPI_Wtime(), t1;
#define LAP(phase) (t1 = MPI_Wtime(), t[phase] += t1 - t0, t0 = t1)

    if (manual) pack_all(d, u);
    LAP(T_PACK);

    if (o->mode == MODE_BLOCKING) {
        /* Shift pattern: send across face f, receive the opposite ghost. */
        for (int f = 0; f < 6; f++) {
            int g = f ^ 1;
            if (manual)
                MPI_Sendrecv(d->sbuf[f], d->fcount[f], MPI_DOUBLE, d->nbr[f], f,
                             d->rbuf[g], d->fcount[g], MPI_DOUBLE, d->nbr[g], f,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            else
                MPI_Sendrecv(u, 1, d->stype[f], d->nbr[f], f,
                             u, 1, d->rtype[g], d->nbr[g], f,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        LAP(T_WAIT);
    } else {
        /* Ghost f comes from nbr[f], which sent it across its face f^1. */
        for (int f = 0; f < 6; f++) {
            if (manual)
                MPI_Irecv(d->rbuf[f], d->fcount[f], MPI_DOUBLE, d->nbr[f], f ^ 1, MPI_COMM_WORLD, &req[f]);
            else
                MPI_Irecv(u, 1, d->rtype[f], d->nbr[f], f ^ 1, MPI_COMM_WORLD, &req[f]);
        }
        for (int f = 0; f < 6; f++) {
            if (manual)
                MPI_Isend(d->sbuf[f], d->fcount[f], MPI_DOUBLE, d->nbr[f], f, MPI_COMM_WORLD, &req[6 + f]);
            else
                MPI_Isend(u, 1, d->stype[f], d->nbr[f], f, MPI_COMM_WORLD, &req[6 + f]);
        }
        LAP(T_POST);
        if (o->mode != MODE_NONBLOCKING) {
            compute_interior(d, u, un, o->r, o->mode == MODE_OVERLAP_TEST ? req : NULL, o->test_every);
            LAP(T_CINT);
        }
        MPI_Waitall(12, req, MPI_STATUSES_IGNORE);
        LAP(T_WAIT);
    }

    if (manual) unpack_all(d, u);
    LAP(T_UNPACK);

    if (o->mode == MODE_OVERLAP || o->mode == MODE_OVERLAP_TEST) {
        compute_boundary(d, u, un, o->r);
        LAP(T_CBD);
    } else {
        stencil_box(d, u, un, o->r, 1, d->n[0], 1, d->n[1], 1, d->n[2]);
        LAP(T_CINT);
    }
#undef LAP
}

static void init_field(const Dom *d, const Opts *o, double *u)
{
    const double h = 1.0 / (o->N + 1);
    for (int k = 1; k <= d->n[2]; k++) {
        double sz = sin(M_PI * h * (d->off[2] + k));
        for (int j = 1; j <= d->n[1]; j++) {
            double syz = sz * sin(M_PI * h * (d->off[1] + j));
            for (int i = 1; i <= d->n[0]; i++)
                u[IDX(d, i, j, k)] = syz * sin(M_PI * h * (d->off[0] + i));
        }
    }
}

/* max |u - g^K u0| / g^K over the whole domain */
static double check_error(const Dom *d, const Opts *o, const double *u, int K)
{
    const double h = 1.0 / (o->N + 1);
    const double s = sin(M_PI * h / 2);
    const double G = pow(1.0 - 12.0 * o->r * s * s, K);
    double err = 0;
    for (int k = 1; k <= d->n[2]; k++) {
        double sz = sin(M_PI * h * (d->off[2] + k));
        for (int j = 1; j <= d->n[1]; j++) {
            double syz = sz * sin(M_PI * h * (d->off[1] + j));
            for (int i = 1; i <= d->n[0]; i++) {
                double e = fabs(u[IDX(d, i, j, k)] - G * syz * sin(M_PI * h * (d->off[0] + i)));
                if (e > err) err = e;
            }
        }
    }
    double gerr;
    MPI_Allreduce(&err, &gerr, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return gerr / G;
}

/* Bytes sent per step by this rank, classified by where the neighbor runs:
   [0] other node, [1] same node / other CCD, [2] same CCD. */
static void placement_stats(const Dom *d, const Opts *o, int *node_of, int *cpu_of,
                            double bytes[3], int *nodes)
{
    MPI_Comm nodecomm;
    int lrank, leader = rank_, cpu = -1;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &nodecomm);
    MPI_Comm_rank(nodecomm, &lrank);
    MPI_Allreduce(MPI_IN_PLACE, &leader, 1, MPI_INT, MPI_MIN, nodecomm);
    MPI_Comm_free(&nodecomm);
#ifdef __linux__
    cpu = sched_getcpu();
#endif
    MPI_Allgather(&leader, 1, MPI_INT, node_of, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&cpu, 1, MPI_INT, cpu_of, 1, MPI_INT, MPI_COMM_WORLD);

    int is_leader = (lrank == 0);
    MPI_Allreduce(&is_leader, nodes, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    bytes[0] = bytes[1] = bytes[2] = 0;
    for (int f = 0; f < 6; f++) {
        int q = d->nbr[f];
        if (q == MPI_PROC_NULL) continue;
        double b = 8.0 * d->fcount[f];
        if (node_of[q] != node_of[rank_]) bytes[0] += b;
        else if (o->ccd_size > 0 && cpu >= 0 && cpu_of[q] >= 0 &&
                 cpu / o->ccd_size == cpu_of[q] / o->ccd_size) bytes[2] += b;
        else bytes[1] += b;
    }
}

static void print_placement(const Dom *d, const int *cpu_of)
{
    char line[160], host[MPI_MAX_PROCESSOR_NAME];
    int len;
    MPI_Get_processor_name(host, &len);
    snprintf(line, sizeof line, "rank %4d  coords (%d,%d,%d)  local %dx%dx%d  host %s  cpu %d",
             rank_, d->c[0], d->c[1], d->c[2], d->n[0], d->n[1], d->n[2], host, cpu_of[rank_]);
    char *all = rank_ == 0 ? malloc(sizeof line * size_) : NULL;
    MPI_Gather(line, sizeof line, MPI_CHAR, all, sizeof line, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (rank_ == 0) {
        for (int p = 0; p < size_; p++) printf("%s\n", all + p * sizeof line);
        free(all);
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &size_);

    Opts o;
    parse(argc, argv, &o);
    if (o.P[0] == 0) {
        /* MPI_Dims_create returns non-increasing dims; put the largest count on z
           and the smallest on x so the unit-stride direction stays long. */
        int dims[3] = {0, 0, 0};
        MPI_Dims_create(size_, 3, dims);
        o.P[0] = dims[2]; o.P[1] = dims[1]; o.P[2] = dims[0];
    }
    if (o.P[0] * o.P[1] * o.P[2] != size_) die("Px*Py*Pz must equal the number of ranks");
    for (int dim = 0; dim < 3; dim++)
        if (o.P[dim] < 1 || o.P[dim] > o.N) die("each Pdim must be between 1 and N");

    Dom d;
    setup_domain(&d, &o);
    double *u = calloc(d.len, sizeof(double)), *un = calloc(d.len, sizeof(double));
    if (!u || !un) die("out of memory (field)");
    init_field(&d, &o, u);

    int *node_of = malloc(sizeof(int) * size_), *cpu_of = malloc(sizeof(int) * size_), nodes;
    double bytes[3];
    placement_stats(&d, &o, node_of, cpu_of, bytes, &nodes);
    if (o.verbose) print_placement(&d, cpu_of);

    double t[T_NUM] = {0};
    for (int s = 0; s < o.warmup; s++) {
        step(&d, &o, u, un, t);
        double *tmp = u; u = un; un = tmp;
    }
    memset(t, 0, sizeof t);

    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();
    for (int s = 0; s < o.steps; s++) {
        step(&d, &o, u, un, t);
        double *tmp = u; u = un; un = tmp;
    }
    double ttot = MPI_Wtime() - t0;

    double rel_err = check_error(&d, &o, u, o.warmup + o.steps);

    /* Per-step phase times: average and max over ranks. */
    double tavg[T_NUM], tmax[T_NUM], tstep;
    for (int p = 0; p < T_NUM; p++) t[p] /= o.steps;
    ttot /= o.steps;
    MPI_Reduce(t, tavg, T_NUM, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(t, tmax, T_NUM, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ttot, &tstep, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    double my_bytes = bytes[0] + bytes[1] + bytes[2], max_bytes, sum_bytes[3];
    int my_msgs = 0, max_msgs;
    for (int f = 0; f < 6; f++) my_msgs += (d.nbr[f] != MPI_PROC_NULL);
    long my_cells = (long)d.n[0] * d.n[1] * d.n[2], max_cells, min_cells;
    MPI_Reduce(&my_bytes, &max_bytes, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(bytes, sum_bytes, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&my_msgs, &max_msgs, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&my_cells, &max_cells, 1, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&my_cells, &min_cells, 1, MPI_LONG, MPI_MIN, 0, MPI_COMM_WORLD);

    const int pass = rel_err < 1e-10;
    if (rank_ == 0) {
        for (int p = 0; p < T_NUM; p++) tavg[p] /= size_;
        printf("halo3d  N=%d  grid=%dx%dx%d  P=%d  nodes=%d  mode=%s  pack=%s  steps=%d\n",
               o.N, o.P[0], o.P[1], o.P[2], size_, nodes, mode_name[o.mode], pack_name[o.pack], o.steps);
        printf("  time/step %.3f ms   (avg | max over ranks, ms)\n", 1e3 * tstep);
        for (int p = 0; p < T_NUM; p++)
            printf("    %-9s %9.4f | %9.4f\n", tname[p], 1e3 * tavg[p], 1e3 * tmax[p]);
        printf("  halo bytes/step: max per rank %.0f; total off-node %.0f, same-node %.0f, same-CCD %.0f\n",
               max_bytes, sum_bytes[0], sum_bytes[1], sum_bytes[2]);
        printf("  local cells min %ld max %ld   rel_err %.3e  %s\n",
               min_cells, max_cells, rel_err, pass ? "PASS" : "FAIL");

        if (o.csv) {
            FILE *fp = fopen(o.csv, "a");
            if (!fp) { perror(o.csv); }
            else {
                fseek(fp, 0, SEEK_END);
                if (ftell(fp) == 0) {
                    fprintf(fp, "tag,tls,jobid,mode,pack,N,Px,Py,Pz,P,nodes,steps,t_step");
                    for (int p = 0; p < T_NUM; p++) fprintf(fp, ",%s_avg,%s_max", tname[p], tname[p]);
                    fprintf(fp, ",halo_bytes_max,msgs_max,bytes_offnode,bytes_samenode,bytes_sameccd,"
                                "cells_min,cells_max,rel_err\n");
                }
                const char *job = getenv("SLURM_JOB_ID");
                /* UCX_TLS is a comma list (sm,self,rc); store it with '+' */
                char tls[128];
                snprintf(tls, sizeof tls, "%s", getenv("UCX_TLS") ? getenv("UCX_TLS") : "-");
                for (char *c = tls; *c; c++) if (*c == ',') *c = '+';
                fprintf(fp, "%s,%s,%s,%s,%s,%d,%d,%d,%d,%d,%d,%d,%.6e", o.tag, tls, job ? job : "-",
                        mode_name[o.mode], pack_name[o.pack], o.N, o.P[0], o.P[1], o.P[2], size_,
                        nodes, o.steps, tstep);
                for (int p = 0; p < T_NUM; p++) fprintf(fp, ",%.6e,%.6e", tavg[p], tmax[p]);
                fprintf(fp, ",%.0f,%d,%.0f,%.0f,%.0f,%ld,%ld,%.3e\n", max_bytes, max_msgs,
                        sum_bytes[0], sum_bytes[1], sum_bytes[2], min_cells, max_cells, rel_err);
                fclose(fp);
            }
        }
    }

    for (int f = 0; f < 6; f++) {
        free(d.sbuf[f]); free(d.rbuf[f]);
        MPI_Type_free(&d.stype[f]); MPI_Type_free(&d.rtype[f]);
    }
    free(u); free(un); free(node_of); free(cpu_of);
    MPI_Finalize();
    return pass ? 0 : 1;
}
