#include <stdio.h>

#include “mpi.h”

int main(int argc, char** argv)

{

  int i, rank, size, dest, src, tag;
  double data[1024];
  MPI_Datatype Data16Type;
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  src = 0;
  dest = 1;
  tag = 0;

  if (!rank) {
    for (i = 0; i < 1024; i++) data[i] = (double)i;
    MPI_Type_vector(64, 1, 16, MPI_DOUBLE, &Data16Type);
    MPI_Type_commit(&Data16Type);
    MPI_Send(data, 1, Data16Type, dest, tag, MPI_COMM_WORLD);
    MPI_Type_free(&Data16Type);
  }  else {
    MPI_Recv(data, 64, MPI_DOUBLE, src, tag, MPI_COMM_WORLD, &status);
    for (i = 0; i < 64; i++) printf(” %lf”, data[i]);
    printf(“\n”);
  }
  MPI_Finalize();
}


/* Result:
0.000000  
 16.000000  
 32.000000 
 48.000000                                        
 64.000000
 80.000000
 96.000000
 112.000000
 ...
*/
