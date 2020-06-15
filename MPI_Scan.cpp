//Source: https://cvw.cac.cornell.edu/MPIcc/scan
// Modified

#include "mpi.h"
#include <math.h>
#define WCOMM MPI_COMM_WORLD
#define NUMPTS 10
main(int argc, char **argv){
  int mype, ierr;
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(WCOMM, &mype);

  double l = 0.5;
  double exp_sum = 0;
  double exp_pdf_i = 0.0;
  double exp_cdf_i = 0.0;
  double DIV_CONST = 2.0;
  int i;

  for(i = 0; i < NUMPTS; i++)
  {
    if (i == mype)
    {
      exp_pdf_i = l*exp(-l * ((double) i) / DIV_CONST);
    }
  }

  ierr = MPI_Scan( \
    &exp_pdf_i, &exp_cdf_i, 1, MPI_DOUBLE, MPI_SUM, WCOMM);

  for (i = 0; i < NUMPTS; i++)
  {
    if (i == mype)
    {
      printf("process %d: cumulative sum = %lf, exp_pdf_i= %lf\n", \
             mype, exp_cdf_i, exp_pdf_i);
    }
  }
  ierr = MPI_Finalize();
}
