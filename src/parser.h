#ifndef PARSER_H
#define PARSER_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef enum bool{
    false=0, 
    true
}bool;

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

typedef struct LogGPO {
    int n_size;
    double latency;
    double os_0;
    double or_0;
    LoggpoPara para[1500];
} LogGPO, *pLogGPO;

#define PROC_N 10
#define BYTE_N 30
#define L1_OP_PARA(op) \
    struct L1Para op[PROC_N]; 
#define L2_OP_PARA(op) \
    struct L2Para op[PROC_N];
typedef struct IMB_BYTES_Para_L1 {
    int bytes;
    double t_avg;
} Para_L1;
typedef struct IMB_BYTES_Para_L2 {
    int bytes;
    double t_min, t_max, t_avg;
} Para_L2;
struct L1Para { 
    int proc; 
    int n_byte; 
    Para_L1 para[BYTE_N]; 
};
struct L2Para { 
    int proc; 
    int n_byte; 
    Para_L2 para[BYTE_N]; 
};
typedef struct IMBPara {
    int n_pingpong;
    L1_OP_PARA(pingpong)
    int n_pingping;
    L1_OP_PARA(pingping)
    int n_sendrecv;
    L2_OP_PARA(sendrecv)
    int n_exchange;
    L2_OP_PARA(exchange)
    int n_allreduce;
    L2_OP_PARA(allreduce)
    int n_reduce;
    L2_OP_PARA(reduce)
    int n_reduce_scatter;
    L2_OP_PARA(reduce_scatter)
    int n_allgather;
    L2_OP_PARA(allgather)
    int n_allgatherv;
    L2_OP_PARA(allgatherv)
    int n_gather;
    L2_OP_PARA(gather)
    int n_gatherv;
    L2_OP_PARA(gatherv)
    int n_scatter;
    L2_OP_PARA(scatter)
    int n_scatterv;
    L2_OP_PARA(scatterv)
    int n_alltoall;
    L2_OP_PARA(alltoall)
    int n_alltoallv;
    L2_OP_PARA(alltoallv)
    int n_bcast;
    L2_OP_PARA(bcast)
    int n_barrier;
    struct Para_L0 {
        int proc; 
        int n_byte;
        double t_min, t_max, t_avg;
    } barrier[PROC_N];
} IMBPara, *pIMBPara;

typedef enum IMBState {
    init = 0,
    pingpong,
    pingping,
    sendrecv,
    exchange,
    allreduce,
    reduce,
    reduce_scatter,
    allgather,
    allgatherv,
    gather,
    gatherv,
    scatter,
    scatterv,
    alltoall,
    alltoallv,
    bcast,
    barrier,
    end,
    data,
    unknown
} IMBState;

void parse_loggpo(const char* f_para, pLogGPO logps);
void parse_imb(const char* f_para, pIMBPara pimb);

#endif
