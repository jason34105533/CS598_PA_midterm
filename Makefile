# On Campus Cluster: make ARCH=-march=znver3
CC     = mpicc
ARCH  ?= -march=native
CFLAGS = -O3 -std=gnu11 -Wall -Wextra $(ARCH)

all: bin/halo3d bin/pingpong bin/halo3d_irr bin/halo3d_ext bin/pingpong_multi

bin/%: src/%.c | bin
	$(CC) $(CFLAGS) -o $@ $< -lm

bin:
	mkdir -p bin

clean:
	rm -rf bin

.PHONY: all clean
