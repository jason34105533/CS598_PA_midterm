/*
 * pingpong.c -- one-way latency and bandwidth between ranks 0 and 1, used to
 * fit the alpha-beta model T(m) = alpha + m * beta for a given placement
 * (same CCD, other CCD, other socket, other node).
 *
 * Each message size is timed in TRIALS separate trials of `reps` round trips;
 * the minimum and the median over trials are reported (the mean of a single
 * trial was too noisy at mid sizes).
 *
 * usage: pingpong [--max-bytes B] [--trials T] [--csv FILE] [--tag STR]
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
    if (size < 2) {
        if (rank == 0) fprintf(stderr, "pingpong needs at least 2 ranks\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    long max_bytes = 1L << 24;
    int trials = 7;
    const char *csv = NULL, *tag = "-";
    for (int a = 1; a + 1 < argc; a += 2) {
        if (!strcmp(argv[a], "--max-bytes")) max_bytes = atol(argv[a + 1]);
        else if (!strcmp(argv[a], "--trials")) trials = atoi(argv[a + 1]);
        else if (!strcmp(argv[a], "--csv")) csv = argv[a + 1];
        else if (!strcmp(argv[a], "--tag")) tag = argv[a + 1];
    }

    char host[MPI_MAX_PROCESSOR_NAME], peer_host[MPI_MAX_PROCESSOR_NAME];
    int len, cpu = -1, peer_cpu = -1;
    MPI_Get_processor_name(host, &len);
#ifdef __linux__
    cpu = sched_getcpu();
#endif
    if (rank == 1) {
        MPI_Send(host, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&cpu, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
    } else if (rank == 0) {
        MPI_Recv(peer_host, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&peer_cpu, 1, MPI_INT, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("pingpong: rank0 %s cpu %d  <->  rank1 %s cpu %d\n", host, cpu, peer_host, peer_cpu);
        const char *tls = getenv("UCX_TLS");
        printf("UCX_TLS=%s  trials=%d\n", tls ? tls : "(unset)", trials);
        printf("%12s %12s %12s %10s %10s\n", "bytes", "min (us)", "median (us)", "GB/s min", "GB/s med");
    }

    char *buf = malloc(max_bytes > 0 ? max_bytes : 1);
    memset(buf, 1, max_bytes > 0 ? max_bytes : 1);
    FILE *fp = NULL;
    if (rank == 0 && csv) {
        fp = fopen(csv, "a");
        if (fp) {
            fseek(fp, 0, SEEK_END);
            if (ftell(fp) == 0)
                fprintf(fp, "tag,tls,jobid,host0,cpu0,host1,cpu1,bytes,reps,trials,"
                            "t_min,t_med,gbps_min,gbps_med\n");
        }
    }

    /* UCX_TLS is a comma list (sm,self,rc); store it with '+' so the CSV stays valid */
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
                if (rank == 0) {
                    MPI_Send(buf, (int)m, MPI_CHAR, 1, 2, MPI_COMM_WORLD);
                    MPI_Recv(buf, (int)m, MPI_CHAR, 1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                } else if (rank == 1) {
                    MPI_Recv(buf, (int)m, MPI_CHAR, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(buf, (int)m, MPI_CHAR, 0, 2, MPI_COMM_WORLD);
                }
            }
            t[tr] = (MPI_Wtime() - t0) / (2.0 * reps);
        }
        if (rank == 0) {
            qsort(t, trials, sizeof(double), cmp_dbl);
            double tmin = t[0];
            double tmed = trials % 2 ? t[trials / 2] : 0.5 * (t[trials / 2 - 1] + t[trials / 2]);
            printf("%12ld %12.3f %12.3f %10.3f %10.3f\n", m, 1e6 * tmin, 1e6 * tmed,
                   m / tmin / 1e9, m / tmed / 1e9);
            if (fp) {
                const char *job = getenv("SLURM_JOB_ID");
                fprintf(fp, "%s,%s,%s,%s,%d,%s,%d,%ld,%d,%d,%.6e,%.6e,%.4f,%.4f\n", tag,
                        tls_csv, job ? job : "-", host, cpu, peer_host, peer_cpu, m,
                        reps, trials, tmin, tmed, m / tmin / 1e9, m / tmed / 1e9);
            }
        }
    }
    free(t);
    if (fp) fclose(fp);
    free(buf);
    MPI_Finalize();
    return 0;
}
