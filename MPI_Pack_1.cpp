float* entries;
int* column_subscripts;
int nonzeroes; /* number of nonzeroes in row */
int position;
int row_number;
char* buffer[HUGE]; /* HUGE is a predefined constant */
MPI_Status status;

if (my_rank == 0) {
  /* Get the number of nonzeros in the row. */
  /* Allocate storage for the row. */
  /* Initialize entries and column_subscripts */

  /* Now pack the data and send */
  position = 0;
  MPI_Pack(&nonzeroes, 1, MPI_INT, buffer, HUGE,
  &position, MPI_COMM_WORLD);
  MPI_Pack(&row_number, 1, MPI_INT, buffer, HUGE,
  &position, MPI_COMM_WORLD);
  MPI_Pack(entries, nonzeroes, MPI_FLOAT, buffer,
  HUGE, &position, MPI_COMM_WORLD);
  MPI_Pack(column_subscripts, nonzeroes, MPI_INT,
  buffer, HUGE, &position, MPI_COMM_WORLD);
  MPI_Send(buffer, position, MPI_PACKED, 1, 193,
  MPI_COMM_WORLD);
} else { /* my_rank == 1 */
  MPI_Recv(buffer, HUGE, MPI_PACKED, 0, 193,
  MPI_COMM_WORLD, &status);
  position = 0;
  MPI_Unpack(buffer, HUGE, &position, &nonzeroes,
  1, MPI_INT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, HUGE, &position, &row_number,
  1, MPI_INT, MPI_COMM_WORLD);
  /* Allocate storage for entries and column_subscripts */
  entries = (float *) malloc(nonzeroes*sizeof(float));
  column_subscripts = (int *) malloc(nonzeroes*sizeof(int));
  MPI_Unpack(buffer,HUGE, &position, entries,
  nonzeroes, MPI_FLOAT, MPI_COMM_WORLD);
  MPI_Unpack(buffer, HUGE, &position, column_subscripts,
  nonzeroes, MPI_INT, MPI_COMM_WORLD);
}
