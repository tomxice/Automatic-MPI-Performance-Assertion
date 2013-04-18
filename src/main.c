#include "mpi.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "timer.h"

#define SIZE (20)
char array[SIZE];

int main( int argc, char** argv) {
    int myrank;
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(array,SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
