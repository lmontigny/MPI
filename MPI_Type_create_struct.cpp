#include <stdio.h>
#include <mpi.h>

int main(int argc, char ** argv) {
    int rank, size;
    char name[80];
    int length;

    MPI_Init(&argc, &argv); 

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Get_processor_name(name,&length);

    typedef struct customData{
        double iv;
        double dv[5];
        char cv[10];
    } customData;

    MPI_Datatype type;
    MPI_Datatype ctype[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_CHAR};
    int blocklen[3] = {1, 5, 10};
    MPI_Aint disp[3];

    disp[0] = offsetof(customData, iv);
    disp[1] = offsetof(customData, dv);
    disp[2] = offsetof(customData, cv);

    MPI_Type_create_struct(3, blocklen, disp, ctype, &type);
    MPI_Type_commit(&type);

    if(rank == 0){
        customData data = {5,
            {1,2,3,4,5}, 
            {'h','d','h','a','q','w','e','s','l','z'}};
        printf("%f %.2f %c\n", data.iv, data.dv[0], data.cv[0]);

        MPI_Send(&data, 1, type, 1, 0, MPI_COMM_WORLD); 
    } else {
        customData recv;
        MPI_Status status;
        MPI_Recv(&recv, 1, type, 0, 0, MPI_COMM_WORLD, &status);
        printf("%f %.2f %c\n", recv.iv, recv.dv[0], recv.cv[0]);
    }

    MPI_Finalize();
}
