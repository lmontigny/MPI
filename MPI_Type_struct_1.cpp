#include <stdio.h>
#include “mpi.h”
#define ROWS  10
#define COLS  10

struct _cell {
  int coord[2];           /* cell coordinates */
  double elevation;       /* cell elevation */
  char landcover;         /* land cover type */
  double albedos[4];      /* seasonal albedos */
};

struct _cell cells[ROWS][COLS];


int main(int argc, char** argv)

{
  int i, j, rank, size, dest, src, tag, base, blocklen[] = {2, 1, 1, 4};

  MPI_Datatype Cells;
  MPI_Datatype type[] = {MPI_INT, MPI_DOUBLE, MPI_CHAR, MPI_DOUBLE};
  MPI_Aint disp[4];
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);


  /* compute displacements of structure components */
  MPI_Address(cells, disp);
  MPI_Address(&cells[0][0].elevation, disp+1);
  MPI_Address(&cells[0][0].landcover, disp+2);
  MPI_Address(&cells[0][0].albedos, disp+3);
  base = disp[0];
  for (i = 0; i < 4; i++) disp[i] -= base;
  MPI_Type_struct(4, blocklen, disp, type, &Cells);

  MPI_Type_commit(&Cells);
  src = 0;
  dest = 1;
  tag = 0;
  if (!rank) {
    for (j = 0; j < ROWS; j++) {
      for (i = 0; i < COLS; i++) {
        cells[j][i].coord[0] = i;
        cells[j][i].coord[1] = j;
        cells[j][i].elevation = (double)(i * 1000 + j);
      }
    }
    printf(“%d: Sending cells to rank %d\n”, rank, dest);
    MPI_Send(cells, (ROWS * COLS), Cells, dest, tag, MPI_COMM_WORLD);
  }  else {
    printf(“%d: Receiving cells from rank %d\n”, rank, src);
    MPI_Recv(cells, (ROWS * COLS), Cells, src, tag, MPI_COMM_WORLD,
&status);
    for (j = 0; j < ROWS; j++)
      for (i = 0; i < COLS; i++)
        printf(“i=%d, j=%d, coords=(%d,%d), elev=%lf\n”,
          i, j, cells[j][i].coord[0],
          cells[j][i].coord[1],
          cells[j][i].elevation);
  }

  MPI_Type_free(&Cells);
  MPI_Finalize();

}
