// ----------------- C++ -----------------

#include <cstring>
#include <mpi.h>

using namespace MPI;

int main(int argc, char **argv) {
    const char *msg="Hello!";
    const int len=strlen(msg)+1;
    char *buf = new char[len];

    Init(argc, argv);
    int rank = COMM_WORLD.Get_rank();

    if (rank == 0) {
        Request isreq = COMM_WORLD.Isend(msg, len, MPI_CHAR, 0, 0);
        Request irreq = COMM_WORLD.Irecv(buf, len, MPI_CHAR, 0, 0);
        isreq.Cancel();
        irreq.Cancel();
    }

    Finalize();
    return 0;
}

// ----------------- C -----------------
#include <stdio.h>
#include <string.h>
#include <mpi.h>

int main(int argc, char **argv) {
    int rank;
    const char *msg="Hello!";
    const int len=strlen(msg)+1;
    char  buf[len];

    MPI_Request isreq, irreq;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        MPI_Isend((void*)msg, len, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &isreq);
        MPI_Irecv(buf, len, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &irreq);
        MPI_Cancel(&irreq);
        MPI_Cancel(&isreq);
    }


    MPI_Finalize();
    return 0;
}
