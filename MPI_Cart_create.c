// Author: Mary Thomas, Department of Computer ScienceComputational Science Research Center (CSRC)San Diego State University (SDSU)

#include <stdio.h>
#include "mpi.h"
#include <math.h>

main(int argc, char *argv[]) {
    int       ndims=1;
    int       p;
    int       my_rank;
    MPI_Comm  comm1D;
    int       dims[ndims],coord[ndims];
    int       wrap_around[ndims];
    int       reorder;
    int       my_cart_rank;
    int       ierr;
    int       nrows, ncols;
    
    /* start up initial MPI environment */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    dims[0]=0;

    /* create cartesian topology for processes */
    MPI_Dims_create(p, ndims, dims);
    if( my_rank == 0 )printf("PW[%d]/[%d%]: PEdims = [%d] \n",my_rank,p,dims[0]);
    
    /* create cartesian mapping */
    wrap_around[0] = 0;
    
    // set periodicity .false.
    reorder = 1;
    ierr =0;
    ierr = MPI_Cart_create(MPI_COMM_WORLD, ndims, dims,wrap_around, reorder, &comm1D);
    if(ierr != 0) printf("ERROR[%d] creating CART\n",ierr);
    
    /* find my coordinates in the cartesian communicator group */
    MPI_Cart_coords(comm1D, my_rank, ndims, coord);
    
    /* use my coordinates to find my rank in cartesian group*/
    MPI_Cart_rank(comm1D, coord, &my_cart_rank);
    printf("PW[%d]: my_cart_rank PCM[%d], my coords = (%d)\n",my_rank, my_cart_rank, coord[0]);
    MPI_Comm_free(&comm1D );
    MPI_Finalize();

}  /* main */
