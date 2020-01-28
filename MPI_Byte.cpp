#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define SIZE 10

struct data_types {
    int a;
    double x;
};

int main(){
    int rank,size;
    MPI_Init(NULL,NULL);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Status status;

    if ( SIZE < size ) exit(1);
   
    struct data_types array_of_structs[SIZE];
    struct data_types single_struct;
   
    if ( 0 == rank )
    for ( int i = 0 ; i < SIZE; i++)
    {
        array_of_structs[i].a = time(NULL);
        array_of_structs[i].x = (i+1)*1.1;
    }

    if ( 0 == rank )
    for ( int i = 0 ; i < size; i++)
    {
        MPI_Send(&(array_of_structs[i]),
            sizeof(struct data_types),MPI_BYTE,i,0,MPI_COMM_WORLD);
    }   

    MPI_Recv(&single_struct,sizeof(struct data_types),MPI_BYTE,0,0,MPI_COMM_WORLD,&status);

    printf("Rank: %d\ta: %d\tx: %f\n",rank,single_struct.a,single_struct.x);
    MPI_Finalize();
    return 0;
}
