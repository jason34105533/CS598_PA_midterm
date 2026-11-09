/*
 * pingpong_multi.c -- point-to-point bandwidth with many pairs active at once.
 * Rank r is paired with r + P/2 (P even), so with block placement the first
 * half of the ranks talks to the second half (other node, other socket, ...).
 * All pairs run at the same time; each trial reports the slowest pair.
 *
 *   --pattern pingpong   one message in flight per pair (as pingpong.c)
 *   --pattern exchange   both ranks of a pair send and receive at once
 *                        (Irecv + Isend + Waitall), like one face of a halo exchange
 *
 * Time and bandwidth are per message, per direction. With P = 2 and pattern
 * pingpong this is the same measurement as pingpong.c.
 *
 * usage: pingpong_multi [--pattern P] [--max-bytes B] [--trials T] [--csv FILE] [--tag STR]
 */
#define _GNU_SOURCE
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef __linux__
#include <sched.h>
#endif

static int cmp_dbl(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size < 2 || size % 2) {
        if (rank == 0) fprintf(stderr, "pingpong_multi needs an even number of ranks\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    long max_bytes = 1L << 24;
    int trials = 7, exchange = 0;
    const char *csv = NULL, *tag = "-";
    for (int a = 1; a + 1 < argc; a += 2) {
        if (!strcmp(argv[a], "--max-bytes")) max_bytes = atol(argv[a + 1]);
        else if (!strcmp(argv[a], "--trials")) trials = atoi(argv[a + 1]);
        else if (!strcmp(argv[a], "--csv")) csv = argv[a + 1];
        else if (!strcmp(argv[a], "--tag")) tag = argv[a + 1];
        else if (!strcmp(argv[a], "--pattern")) exchange = !strcmp(argv[a + 1], "exchange");
    }

    const int half = size / 2, pairs = half;
    const int partner = rank < half ? rank + half : rank - half;
    const int leader = rank < half;

    char host[MPI_MAX_PROCESSOR_NAME], peer_host[MPI_MAX_PROCESSOR_NAME];
    int len, cpu = -1, peer_cpu = -1;
    MPI_Get_processor_name(host, &len);
#ifdef __linux__
    cpu = sched_getcpu();
#endif
    if (rank == half) {
        MPI_Send(host, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&cpu, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
    } else if (rank == 0) {
        MPI_Recv(peer_host, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, half, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&peer_cpu, 1, MPI_INT, half, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("pingpong_multi: %d pairs, %s; pair 0: %s cpu %d <-> %s cpu %d\n", pairs,
               exchange ? "exchange" : "pingpong", host, cpu, peer_host, peer_cpu);
        const char *tls = getenv("UCX_TLS");
        printf("UCX_TLS=%s  trials=%d  (times: slowest pair)\n", tls ? tls : "(unset)", trials);
        printf("%12s %12s %12s %10s %10s\n", "bytes", "min (us)", "median (us)", "GB/s min", "GB/s med");
    }

    size_t nb = max_bytes > 0 ? (size_t)max_bytes : 1;
    char *sbuf = malloc(nb), *rbuf = malloc(nb);
    memset(sbuf, 1, nb);
    memset(rbuf, 0, nb);
    FILE *fp = NULL;
    if (rank == 0 && csv) {
        fp = fopen(csv, "a");
        if (fp) {
            fseek(fp, 0, SEEK_END);
            if (ftell(fp) == 0)
                fprintf(fp, "tag,tls,jobid,pattern,pairs,host0,cpu0,host1,cpu1,bytes,reps,trials,"
                            "t_min,t_med,gbps_min,gbps_med\n");
        }
    }

    char tls_csv[128];
    snprintf(tls_csv, sizeof tls_csv, "%s", getenv("UCX_TLS") ? getenv("UCX_TLS") : "-");
    for (char *c = tls_csv; *c; c++) if (*c == ',') *c = '+';

    double *t = malloc(sizeof(double) * (trials > 0 ? trials : 1));
    for (long m = 8; m <= max_bytes; m *= 2) {
        int reps = (int)(((1L << 28) / m) > 2000 ? 2000 : ((1L << 28) / m));
        if (reps < 20) reps = 20;
        for (int tr = 0; tr < trials; tr++) {
            double t0 = 0;
            for (int it = -10; it < reps; it++) { /* 10 warmup iterations per trial */
                if (it == 0) { MPI_Barrier(MPI_COMM_WORLD); t0 = MPI_Wtime(); }
                if (exchange) {
                    MPI_Request rq[2];
                    MPI_Irecv(rbuf, (int)m, MPI_CHAR, partner, 2, MPI_COMM_WORLD, &rq[0]);
                    MPI_Isend(sbuf, (int)m, MPI_CHAR, partner, 2, MPI_COMM_WORLD, &rq[1]);
                    MPI_Waitall(2, rq, MPI_STATUSES_IGNORE);
                } else if (leader) {
                    MPI_Send(sbuf, (int)m, MPI_CHAR, partner, 2, MPI_COMM_WORLD);
                    MPI_Recv(rbuf, (int)m, MPI_CHAR, partner, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                } else {
                    MPI_Recv(rbuf, (int)m, MPI_CHAR, partner, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(sbuf, (int)m, MPI_CHAR, partner, 2, MPI_COMM_WORLD);
                }
            }
            double mine = (MPI_Wtime() - t0) / (exchange ? reps : 2.0 * reps);
            MPI_Allreduce(&mine, &t[tr], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        }
        if (rank == 0) {
            qsort(t, trials, sizeof(double), cmp_dbl);
            double tmin = t[0];
            double tmed = trials % 2 ? t[trials / 2] : 0.5 * (t[trials / 2 - 1] + t[trials / 2]);
            printf("%12ld %12.3f %12.3f %10.3f %10.3f\n", m, 1e6 * tmin, 1e6 * tmed,
                   m / tmin / 1e9, m / tmed / 1e9);
            if (fp) {
                const char *job = getenv("SLURM_JOB_ID");
                fprintf(fp, "%s,%s,%s,%s,%d,%s,%d,%s,%d,%ld,%d,%d,%.6e,%.6e,%.4f,%.4f\n", tag,
                        tls_csv, job ? job : "-", exchange ? "exchange" : "pingpong", pairs, host,
                        cpu, peer_host, peer_cpu, m, reps, trials, tmin, tmed, m / tmin / 1e9,
                        m / tmed / 1e9);
            }
        }
    }
    free(t);
    if (fp) fclose(fp);
    free(sbuf);
    free(rbuf);
    MPI_Finalize();
    return 0;
}
