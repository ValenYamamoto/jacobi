jacobi_exe: pjacobi.c
	gcc -Wall -O2 pjacobi.c matrixUtil.c -lm -pthread

run:
	numactl -C 0-31 ./a.out 32
