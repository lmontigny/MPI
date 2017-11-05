 // Layout of Type_struct
 
 struct my_struct {
    int int_values[10];
    double average;
    char debug_name[MAX_NAME_LEN];
    int flag;
 };
 
void make_datatype(MPI_Datatype *new_type) {
    struct my_struct foo;
    int i, counts[4] = { 10, 1, MAX_NAME_LEN, 1 };
   MPI_Datatype types[4] = { MPI_INT, MPI_DOUBLE, MPI_CHAR, MPI_INT };
   MPI_Aint disps[4];
   MPI_Address(foo.int_values, &disps[0]);
   MPI_Address(&foo.average, &disps[1]);
   MPI_Address(foo.debug_name, &disps[2]);
   MPI_Address(&foo.flag, &disps[3]);
   for (i = 3; i >= 0; --i)
     disps[i] -= disps[0];
   MPI_Type_struct(4, counts, disps, types, new_type);
   MPI_Type_commit(new_type);
}
