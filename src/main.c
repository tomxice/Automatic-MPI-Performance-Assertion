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
    int a = 0;
    int LOOP = 100;
    for (int k = 0; k < LOOP; ++k) {
        //work load
        int W = 10000000;
        for (int i = 0; i < W; ++ i) {
            a = a+i;
        }
        //end work load
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Finalize();
    printf("a:%d\n",a);
    return 0;
}
