#include <iomanip>
#include <iostream>
#include <mpi.h>

using namespace std;

int main(int argc, char** argv)
{
    // Building a Buffer Using MPI_Pack
    struct my_struct {
        int int_value[10];
        double average;
        char debug_name[MAX_NAME_LEN];
        int flag;
    };  
    void send_data(struct my_struct data, MPI_Comm comm, int rank) {
        char buf[BUFSIZE];
        int pos = 0;
        
        MPI_Pack(&data.int_value, 10, MPI_INT, buf, BUFSIZE, &pos, comm);
        MPI_Pack(&data.average, 1, MPI_DOUBLE, buf, BUFSIZE, &pos, comm);
        MPI_Pack(&data.debug_name,MAX_NAME_LEN, MPI_CHAR, buf, BUFSIZE, &pos, comm);
        MPI_Pack(&data.flag, 1, MPI_INT, buf, BUFSIZE, &pos, comm);
        
        MPI_Send(buf, pos, MPI_PACKED, rank, 0, comm);
    }   
 }  
