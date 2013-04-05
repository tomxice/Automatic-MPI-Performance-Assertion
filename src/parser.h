#ifndef PARSER_H
#define PARSER_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef enum bool{
    false=0, 
    true
}bool;

#define VERIFY 0

#define MPS 1500 // Max step in loggpo
typedef struct LoggpoPara {
    int size;
    double os;
    double or;
    double ov;
    double sr;
    double gap;
    double rtt;
    double rtt100;
} LoggpoPara, *pLoggpoPara;

#define PROC_N 10
#define BYTE_N 30
typedef struct IMB_BYTES_Para {
    int bytes;
    double t_min, t_max, t_avg;
} Para;
typedef struct IMBPara {
    int n_barrier;
    int n_bcast;
    int n_reduce;
    int n_gather;
    int n_allreduce;
    int n_allgather;
    struct BarrierPara {
        int proc; 
        int n_byte;
        double t_min, t_max, t_avg;
    } barrier[PROC_N];
    struct BcastPara {
        int proc;
        int n_byte;
        Para para[BYTE_N];
    } bcast[PROC_N];
    struct ReducePara {
        int proc;
        int n_byte;
        Para para[BYTE_N];
    } reduce[PROC_N];
    struct GatherPara {
        int proc;
        int n_byte;
        Para para[BYTE_N];
    } gather[PROC_N];
    struct AllreducePara {
        int proc;
        int n_byte;
        Para para[BYTE_N];
    } allreduce[PROC_N];
    struct AllgatherPara {
        int proc;
        int n_byte;
        Para para[BYTE_N];
    } allgather[PROC_N];
} IMBPara, *pIMBPara;

typedef enum IMBState {
    init,
    barrier,
    bcast,
    reduce,
    gather,
    allreduce,
    allgather,
    end,
    data,
    unknown
} IMBState;

void parse_loggpo(const char* f_para, pLoggpoPara logps);
void parse_imb(const char* f_para, pIMBPara pimb);

#endif
