run.x: milk_crown.c
	qcc -source -D_MPI=1 -DTRACE=2 milk_crown.c
	mpicc -Wall -O2 -std=c99 _milk_crown.c -o run.x -lm

