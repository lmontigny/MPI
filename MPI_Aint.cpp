/**
 * @author RookieHPC
 * @brief Original source code at https://www.rookiehpc.com/mpi/docs/mpi_aint.php
 **/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
 
/**
 * @brief Illustrate how to manipulate the MPI_Aint datatype.
 * @details This application takes the address of elements at different
 * locations and calculates the distance, in bytes, between the two.
 **/
int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
 
    int a[10];
 
    MPI_Aint addr_1;
    MPI_Get_address(&a[2], &addr_1);
 
    MPI_Aint addr_2;
    MPI_Get_address(&a[8], &addr_2);
 
    MPI_Aint addr_gap;
    addr_gap = MPI_Aint_diff(addr_2, addr_1);
 
    printf("Difference between the address of the 3rd int and 9th int is %ld bytes.\n", addr_gap);
 
    MPI_Finalize();
 
    return EXIT_SUCCESS;
}

