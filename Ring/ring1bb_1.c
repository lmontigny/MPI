// Source: perdacherMartin githubgist

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define MASTER 0
#define REPETITIONS 10
#define T_mpisendtype MPI_CHAR

typedef unsigned char T_arraytype;

void InitArray(T_arraytype array[], long size, T_arraytype value );
void Sum(T_arraytype arrayA[], T_arraytype arrayB[], long size);

int main(int argc, char *argv[]){
	int myrank, nprocs, leftid, rightid, val, tmp,i;
	long n;
	double times[REPETITIONS];
	double time,slowest;
	
  	MPI_Status recv_status;
  	MPI_Init(&argc, &argv);
  	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);  
  	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  
	if ( argc != 2 ){
    	printf("Usage: mpiexec -n nodecount ring1bb n \n\n");
		return 1;
    }

	n = atol(argv[1]);

	T_arraytype a[n];
	T_arraytype b[n];
	T_arraytype recv[n];
	InitArray(a,n,myrank);
	
	for ( i = 0 ; i < REPETITIONS ; ++i ){
		MPI_Barrier(MPI_COMM_WORLD);
		
		InitArray(b,n,0);
		
		time = MPI_Wtime();
		
		leftid = myrank - 1 < 0 ? nprocs - 1 : myrank -1 ;
		rightid = myrank + 1 >= nprocs ? 0 : myrank +1 ;

		// send to the left and receive from the righ
		MPI_Sendrecv(&a,n,T_mpisendtype,leftid,0,&recv,n,MPI_INT,rightid,0,MPI_COMM_WORLD, &recv_status);
		Sum(recv, b, n);

		// send to the right and receive from the left
		MPI_Sendrecv(&a,n,T_mpisendtype,rightid,0,&recv,n,MPI_INT,leftid,0,MPI_COMM_WORLD, &recv_status);
		Sum(recv, b, n);
		Sum(a, b, n);
		
		time = MPI_Wtime() - time;
		MPI_Barrier(MPI_COMM_WORLD); 
		// take the slowest time of all processes
		MPI_Reduce(&time, &slowest, 1, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
		
		if (myrank == MASTER) {
	    	times[i]=slowest;
	   	}
	}
	
	if ( myrank == MASTER ){
		double min = times[0];
		// take the best run of the REPITITIONS-runs
	   	for (int i = 1; i < REPETITIONS ; i++) {
			min = ( times[i] < min ) ? times[i] : min;
	   	}
	   	
		printf("%ld:%f\n", n, min);
	}

  	MPI_Finalize();		
}


void InitArray(T_arraytype array[], long size, T_arraytype value ){
	long i;
	
	// initialisation of the array
	for ( i=0l ; i < size ; ++i ){
		array[i] = value;
	}
}

// elementwise sum of A + B and store it in B
void Sum(T_arraytype arrayA[], T_arraytype arrayB[], long size){
	long i;
	
	// calculate elmentwise sum
	for ( i=0l ; i < size ; ++i ){
		arrayB[i] = arrayB[i] + arrayA[i];
	}
}
