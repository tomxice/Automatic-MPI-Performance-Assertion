/********************************************************************
*                                                                   *  
* This file use to estimate how long will a E_MPI call take?
*                                                                   *
********************************************************************/
#include <mpi.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include "estimator.h"
#include "parser.h"


LogGPO log_cmp, log_net, log_smp;
IMBPara imb;
pLocation location;

/********************************************************************
* Determining whether two MPI rank are in the same node or same CPU
********************************************************************/
#define NETWORK 0 //net_para
#define INTRA_NODE 1 //smp_para
#define INTRA_CPU 2 //cmp_para

// TODO you need to modify this function to fit your own arch
int diff_cpu( int c0, int c1 ) {
    if (c0/6 == c1/6) return 0;
    else return 1;
}

int rank_dist( int r0, int r1 ) {
    if (location[r0].node != location[r1].node) return NETWORK;
    else if (diff_cpu(location[r0].core, location[r1].core)) return INTRA_NODE;
    else return INTRA_CPU;
}

pLogGPO getGPO(int r0, int r1) {
    int dist = rank_dist(r0,r1);
    if (dist == NETWORK) return &log_net;
    else if (dist == INTRA_NODE) return &log_smp;
    else if (dist == INTRA_CPU) return &log_cmp;
    else return NULL;
}

/********************************************************************
* For Caculating Time
********************************************************************/
int logbs(pLoggpoPara A, int key, int len) {
    int imin = 0, imax = len-1;
    if (imax < 0) return -2;
    if (key < A[imin].size || key > A[imax].size) return -1;
    // continually narrow search until just one element remains
    while (imin < imax)
    {
        int imid = (imin+imax) >> 1;
        if (imid >= imax) {
            printf ("imid >= imax!\n");
        }
        if (A[imid].size < key)
            imin = imid + 1;
        else
            imax = imid;
    }
    if (A[imin].size > key) --imin;
    return imin;
}

LoggpoPara getPara(pLogGPO p, int size) {
    int a = logbs(p->para, size, p->n_size);
    int b = 0;
    if (a == -1) { // size too large
        a = p->n_size - p->n_size/10;
        b = p->n_size - 1;
    }
    else if (p->para[a].size == size) // hit
        return p->para[a];
    else { // interval
        b = a + 1;
    }
    LoggpoPara np;
    double delta = ((double)size-p->para[a].size)/(p->para[b].size-p->para[a].size);
    np.size = size;
    np.os = p->para[a].os + delta*(p->para[b].os - p->para[a].os);
    np.or = p->para[a].or + delta*(p->para[b].or - p->para[a].or);
    np.ov = p->para[a].ov + delta*(p->para[b].ov - p->para[a].ov);
    np.sr = p->para[a].sr + delta*(p->para[b].sr - p->para[a].sr);
    np.gap = p->para[a].gap + delta*(p->para[b].gap - p->para[a].gap);
    np.rtt = p->para[a].rtt + delta*(p->para[b].rtt - p->para[a].rtt);
    np.rtt100 = p->para[a].rtt100 + delta*(p->para[b].rtt100 - p->para[a].rtt100);
    return np;
}


int E_count2byte ( MPI_Datatype datatype, int count ) {
    int nb = 0;
    PMPI_Type_size(datatype, &nb);
    nb *= count;
    return nb;
}
#define E_GET_COLL_AVG_TIME(n_coll_op, coll_op) \
    int pdown=0, pup=0, bdown=0, bup=0;\
    int hit = 0; \
    for (pdown = 0; pdown < imb.n_coll_op-1; ++ pdown) {\
        pup = pdown + 1; \
        if (imb.coll_op[pdown].proc == np) { \
            pup = pdown; \
            hit = 1; \
            break; \
        } \
        else if (imb.coll_op[pdown].proc < np && np < imb.coll_op[pup].proc) { \
            hit = 1; \
            break; \
        } \
        else if (imb.coll_op[pup].proc == np) { \
            pdown = pup; \
            hit = 1; \
            break; \
        } \
    } \
    if (hit == 0) { \
        return -1; \
    } \
    hit = 0; \
    for (bdown = 0; bdown < imb.coll_op[pdown].n_byte; ++ bdown) {\
        bup = bdown + 1; \
        if (imb.coll_op[pdown].para[bdown].bytes == nb) {\
            bup = bdown;\
            hit = 1;\
            break;\
        }\
        else if (imb.coll_op[pdown].para[bdown].bytes < nb && nb < imb.coll_op[pup].para[bup].bytes) {\
            hit = 1;\
            break;\
        }\
        else if (imb.coll_op[pup].para[bup].bytes == nb) {\
            bdown = bup;\
            hit = 1;\
            break;\
        }\
    }\
    if (hit == 0) { \
        bup = imb.coll_op[pdown].n_byte - 1; \
        bdown = bup - 1; \
    } \
    double delta_t, delta_x, slope;\
    double t1[4] = {imb.coll_op[pdown].para[bdown].t_avg, imb.coll_op[pdown].para[bdown].t_avg,\
                    imb.coll_op[pup].para[bdown].t_avg,   imb.coll_op[pdown].para[bup].t_avg};\
    double t2[4] = {imb.coll_op[pup].para[bdown].t_avg,   imb.coll_op[pdown].para[bup].t_avg,\
                    imb.coll_op[pup].para[bup].t_avg,     imb.coll_op[pup].para[bup].t_avg};\
    double x1[4] = {imb.coll_op[pdown].proc,              imb.coll_op[pdown].para[bdown].bytes,\
                    imb.coll_op[pup].para[bdown].bytes,   imb.coll_op[pdown].proc};\
    double x2[4] = {imb.coll_op[pup].proc,                imb.coll_op[pdown].para[bup].bytes,\
                    imb.coll_op[pup].para[bup].bytes,     imb.coll_op[pup].proc};\
    double x[4] = { np, nb, nb, np };\
    double retVal = 0;\
    if (pup == pdown && bup == bdown) { \
        retVal = t1[0]; \
    } \
    else if (pup == pdown) { \
        delta_t = t1[1] - t2[1]; \
        delta_x = x1[1] - x2[1]; \
        slope = delta_t/delta_x; \
        retVal = t1[1] + slope*(x[1]-x1[1]); \
    } \
    else if (bup == bdown) { \
        delta_t = t1[0] - t2[0]; \
        delta_x = x1[0] - x2[0]; \
        slope = delta_t/delta_x; \
        retVal = t1[0] + slope*(x[0]-x1[0]); \
    } \
    else { \
        for (int i = 0; i < 4; ++ i) { \
            delta_t = t1[i] - t2[i]; \
            delta_x = x1[i] - x2[i]; \
            slope = delta_t/delta_x; \
            retVal += t1[i] + slope*(x[i]-x1[i]); \
        } \
        retVal /= 4; \
    }\
    return retVal;

/******************************************************************
*                                                                 *
*                   MPI Functions for Management                  *
*                                                                 *
******************************************************************/
double E_MPI_Init(int * argc, char*** argv)
{
    // assume all data files are existing
    // users may run IMB manually.
    // and copy datas to all machines manually
    parse_loggpo("paras/cmp_para", &log_cmp);
    parse_loggpo("paras/net_para", &log_net);
    parse_loggpo("paras/smp_para", &log_smp);
    parse_imb("paras/coll_para", &imb);
    // get self location HOSTNAME:CORE
    char proc_file_name[50];
    sprintf(proc_file_name, "/proc/%d/stat", getpid());
    FILE* proc_file = fopen(proc_file_name, "r");
    if (! proc_file) {
        printf("Proc File %s Open Failed!\n", proc_file_name);
    }
    int core;
    if (1 != fscanf(proc_file, "%*d %*s %*c %*d %*d %*d %*d %*d %*u %*u %*u %*u %*u %*u %*u %*d %*d %*d %*d %*d %*d %*u %*u %*d %*u %*u %*u %*u %*u %*u %*u %*u %*u %*u %*u %*u %*u %*d %d %*u %*u %*u",&core))
        printf("Read Core ID Failed!\n");
    char hostname[40];
    gethostname(hostname,40);
    // Send their info to rank:0
    int gsize;
    PMPI_Comm_size( MPI_COMM_WORLD, &gsize);
    location = (pLocation)malloc(gsize*sizeof(Location));
    int myrank;
    PMPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank != 0) {
        char sendbuf[100];
        sprintf(sendbuf, "%s %d",hostname,core);
        PMPI_Send( sendbuf, 100, MPI_CHAR, 0, myrank, MPI_COMM_WORLD);
    }
    if (myrank == 0) {
        char rbuf[100];
        MPI_Status ms;
        char** ls = (char**)malloc(gsize*sizeof(char*));
        int cnt = 0;
        ls[0] = (char*)malloc(40);
        if (ls[0] == strcpy(ls[0],hostname)) 
            ++ cnt;
        location[0].node = 0;
        location[0].core = core;
        char r_hn[40];
        int r_core;
        for (int rank = 1; rank < gsize; ++ rank) {
            PMPI_Recv(rbuf, 100, MPI_CHAR, rank, rank, MPI_COMM_WORLD, &ms);
            sscanf(rbuf,"%s %d",r_hn,&r_core);
            int i;
            for (i = 0; i < cnt; ++ i) {
                if(strcmp(ls[i],r_hn) == 0) {
                    location[rank].node = i;
                    location[rank].core = r_core;
                    break;
                }
            }
            if (i == cnt) {
                ls[i] = (char*)malloc(40);
                if (ls[i] == strcpy(ls[i],r_hn)) 
                    ++ cnt;
                location[rank].node = i;
                location[rank].core = r_core;
            }
        }
#if 0
        for (int rank = 0; rank < gsize; ++ rank) {
            printf("rank:%d in node:%d core:%d\n", rank, location[rank].node, location[rank].core);
        }
#endif
    }
#if 0
    LoggpoPara para = getPara(getGPO(0,1),2<<20);
    pLoggpoPara p = &para; 
    printf("\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
            p->size, p->os, p->or, p->ov, p->sr, p->gap, p->rtt, p->rtt100);
#endif
	return 0;
}
#ifdef PERF_MPI_THREADED
double E_MPI_Init_thread (argc, argv, required, provided )
int * argc;
char *** argv;
int required;
int *provided;
{
	return 0;
}
#endif /* PERF_MPI_THREADED */
double E_MPI_Finalize(  )
{
	return 0;
}
/******************************************************************
*                                                                 *
*          MPI Point-to-Point(P2P) Communication Related          *
*                                                                 *
******************************************************************/
double E_MPI_Bsend( buf, count, datatype, dest, tag, comm )
void * buf;
int count;
MPI_Datatype datatype;
int dest;
int tag;
MPI_Comm comm;
{
	return 0;
}
double E_MPI_Bsend_init( buf, count, datatype, dest, tag, comm, request )
void * buf;
int count;
MPI_Datatype datatype;
int dest;
int tag;
MPI_Comm comm;
MPI_Request * request;
{
	return 0;
}
double E_MPI_Recv_init( buf, count, datatype, source, tag, comm, request )
void * buf;
int count;
MPI_Datatype datatype;
int source;
int tag;
MPI_Comm comm;
MPI_Request * request;
{
	return 0;
}
double E_MPI_Send_init( buf, count, datatype, dest, tag, comm, request )
void * buf;
int count;
MPI_Datatype datatype;
int dest;
int tag;
MPI_Comm comm;
MPI_Request * request;
{
	return 0;
}
double E_MPI_Ibsend( buf, count, datatype, dest, tag, comm, request )
void * buf;
int count;
MPI_Datatype datatype;
int dest;
int tag;
MPI_Comm comm;
MPI_Request * request;
{
	return 0;
}
double E_MPI_Irecv( buf, count, datatype, source, tag, comm, request )
void * buf;
int count;
MPI_Datatype datatype;
int source;
int tag;
MPI_Comm comm;
MPI_Request * request;
{
	return 0;
}
double E_MPI_Irsend( buf, count, datatype, dest, tag, comm, request )
void * buf;
int count;
MPI_Datatype datatype;
int dest;
int tag;
MPI_Comm comm;
MPI_Request * request;
{
	return 0;
}
double E_MPI_Isend( buf, count, datatype, dest, tag, comm, request )
void * buf;
int count;
MPI_Datatype datatype;
int dest;
int tag;
MPI_Comm comm;
MPI_Request * request;
{
	return 0;
}
double E_MPI_Issend( buf, count, datatype, dest, tag, comm, request )
void * buf;
int count;
MPI_Datatype datatype;
int dest;
int tag;
MPI_Comm comm;
MPI_Request * request;
{
	return 0;
}
double E_MPI_Recv( buf, count, datatype, source, tag, comm, status )
void * buf;
int count;
MPI_Datatype datatype;
int source;
int tag;
MPI_Comm comm;
MPI_Status * status;
{
	return 0;
}
double E_MPI_Rsend( buf, count, datatype, dest, tag, comm )
void * buf;
int count;
MPI_Datatype datatype;
int dest;
int tag;
MPI_Comm comm;
{
	return 0;
}
double E_MPI_Rsend_init( buf, count, datatype, dest, tag, comm, request )
void * buf;
int count;
MPI_Datatype datatype;
int dest;
int tag;
MPI_Comm comm;
MPI_Request * request;
{
	return 0;
}
double E_MPI_Send( buf, count, datatype, dest, tag, comm )
void * buf;
int count;
MPI_Datatype datatype;
int dest;
int tag;
MPI_Comm comm;
{
	return 0;
}
double E_MPI_Sendrecv( sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status )
void * sendbuf;
int sendcount;
MPI_Datatype sendtype;
int dest;
int sendtag;
void * recvbuf;
int recvcount;
MPI_Datatype recvtype;
int source;
int recvtag;
MPI_Comm comm;
MPI_Status * status;
{
	return 0;
}
double E_MPI_Sendrecv_replace( buf, count, datatype, dest, sendtag, source, recvtag, comm, status )
void * buf;
int count;
MPI_Datatype datatype;
int dest;
int sendtag;
int source;
int recvtag;
MPI_Comm comm;
MPI_Status * status;
{
	return 0;
}
double E_MPI_Ssend( buf, count, datatype, dest, tag, comm )
void * buf;
int count;
MPI_Datatype datatype;
int dest;
int tag;
MPI_Comm comm;
{
	return 0;
}
double E_MPI_Ssend_init( buf, count, datatype, dest, tag, comm, request )
void * buf;
int count;
MPI_Datatype datatype;
int dest;
int tag;
MPI_Comm comm;
MPI_Request * request;
{
	return 0;
}
double E_MPI_Test( request, flag, status )
MPI_Request * request;
int * flag;
MPI_Status * status;
{
	return 0;
}
double E_MPI_Testall( count, array_of_requests, flag, array_of_statuses )
int count;
MPI_Request * array_of_requests;
int * flag;
MPI_Status * array_of_statuses;
{
	return 0;
}
double E_MPI_Testany( count, array_of_requests, index, flag, status )
int count;
MPI_Request * array_of_requests;
int * index;
int * flag;
MPI_Status * status;
{
	return 0;
}
double E_MPI_Test_cancelled( status, flag )
MPI_Status * status;
int * flag;
{
	return 0;
}
double E_MPI_Testsome( incount, array_of_requests, outcount, array_of_indices, array_of_statuses )
int incount;
MPI_Request * array_of_requests;
int * outcount;
int * array_of_indices;
MPI_Status * array_of_statuses;
{
	return 0;
}
double E_MPI_Wait( request, status )
MPI_Request * request;
MPI_Status * status;
{
	return 0;
}
double E_MPI_Waitall( count, array_of_requests, array_of_statuses )
int count;
MPI_Request * array_of_requests;
MPI_Status * array_of_statuses;
{
	return 0;
}
double E_MPI_Waitany( count, array_of_requests, index, status )
int count;
MPI_Request * array_of_requests;
int * index;
MPI_Status * status;
{
	return 0;
}
double E_MPI_Waitsome( incount, array_of_requests, outcount, array_of_indices, array_of_statuses )
int incount;
MPI_Request * array_of_requests;
int * outcount;
int * array_of_indices;
MPI_Status * array_of_statuses;
{
	return 0;
}
double E_MPI_Cancel( request )
MPI_Request * request;
{
	return 0;
}
double E_MPI_Request_free( request )
MPI_Request * request;
{
	return 0;
}
double E_MPI_Start( request )
MPI_Request * request;
{
	return 0;
}
double E_MPI_Startall( count, array_of_requests )
int count;
MPI_Request * array_of_requests;
{
	return 0;
}
double E_MPI_Iprobe( source, tag, comm, flag, status )
int source;
int tag;
MPI_Comm comm;
int * flag;
MPI_Status * status;
{
	return 0;
}
double E_MPI_Probe( source, tag, comm, status )
int source;
int tag;
MPI_Comm comm;
MPI_Status * status;
{
	return 0;
}
/******************************************************************
*                                                                 *
*               MPI Collective Communication Related              *
*                                                                 *
******************************************************************/
double E_MPI_Allgather( sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm )
void * sendbuf;
int sendcount;
MPI_Datatype sendtype;
void * recvbuf;
int recvcount;
MPI_Datatype recvtype;
MPI_Comm comm;
{
    int np, nb; 
    PMPI_Comm_size( comm, &np);
    nb = E_count2byte( sendtype, sendcount);
    E_GET_COLL_AVG_TIME(n_allgather,allgather);
}
double E_MPI_Allgatherv( sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm )
void * sendbuf;
int sendcount;
MPI_Datatype sendtype;
void * recvbuf;
int * recvcounts;
int * displs;
MPI_Datatype recvtype;
MPI_Comm comm;
{
    int np, nb; 
    PMPI_Comm_size( comm, &np);
    nb = E_count2byte( sendtype, sendcount);
    E_GET_COLL_AVG_TIME(n_allgatherv,allgatherv);
}
double E_MPI_Allreduce( sendbuf, recvbuf, count, datatype, op, comm )
void * sendbuf;
void * recvbuf;
int count;
MPI_Datatype datatype;
MPI_Op op;
MPI_Comm comm;
{
    int np, nb; 
    PMPI_Comm_size( comm, &np);
    nb = E_count2byte( datatype, count);
    E_GET_COLL_AVG_TIME(n_allreduce,allreduce);
}
double E_MPI_Alltoall( sendbuf, sendcount, sendtype, recvbuf, recvcnt, recvtype, comm )
void * sendbuf;
int sendcount;
MPI_Datatype sendtype;
void * recvbuf;
int recvcnt;
MPI_Datatype recvtype;
MPI_Comm comm;
{
    int np, nb; 
    PMPI_Comm_size( comm, &np);
    nb = E_count2byte( sendtype, sendcount);
    E_GET_COLL_AVG_TIME(n_alltoall,alltoall);
}
double E_MPI_Alltoallv( sendbuf, sendcnts, sdispls, sendtype, recvbuf, recvcnts, rdispls, recvtype, comm )
void * sendbuf;
int * sendcnts;
int * sdispls;
MPI_Datatype sendtype;
void * recvbuf;
int * recvcnts;
int * rdispls;
MPI_Datatype recvtype;
MPI_Comm comm;
{
    int np, nb; 
    PMPI_Comm_size( comm, &np);
    nb = E_count2byte( sendtype, *sendcnts);
    E_GET_COLL_AVG_TIME(n_alltoallv,alltoallv);
}
double E_MPI_Barrier( comm )
MPI_Comm comm;
{
    int np;
    PMPI_Comm_size( comm, &np);
    int down,up;
    int hit = 0;
    //printf ("n_barrier:%d\n", imb.n_barrier);
    for (down = 0; down < imb.n_barrier-1; ++ down) {
        up = down + 1;
        //printf("np:%d,down:%d,pdown:%d,up:%d,pup:%d\n",np,down,imb.barrier[down].proc,up,imb.barrier[up].proc);
        if (np == imb.barrier[down].proc) {
            up = down;
            hit = 1;
            //printf("hit lower bound\n");
            break;
        }
        else if (imb.barrier[down].proc < np && np < imb.barrier[up].proc) {
            hit = 1;
            //printf("hit internal\n");
            break;
        }
        else if (np == imb.barrier[up].proc) {
            down = up;
            hit = 1;
            //printf("hit upper bound\n");
            break;
        }
    }
    // too many procs
    if (hit == 0) {
        return -1;
    }
    // normal
    double retVal = 0;
    if (down == up) {
        retVal = imb.barrier[up].t_avg;
    }
    else {
        double delta_t = imb.barrier[up].t_avg - imb.barrier[down].t_avg;
        double delta_p = imb.barrier[up].proc - imb.barrier[down].proc;
        double slope = delta_t/delta_p;
        retVal = imb.barrier[down].t_avg + slope*(np - imb.barrier[down].proc);
    }
	return retVal;
}
double E_MPI_Bcast( buffer, count, datatype, root, comm )
void * buffer;
int count;
MPI_Datatype datatype;
int root;
MPI_Comm comm;
{
    int np, nb; 
    PMPI_Comm_size( comm, &np);
    nb = E_count2byte( datatype, count);
    E_GET_COLL_AVG_TIME(n_bcast,bcast)
}
double E_MPI_Gather( sendbuf, sendcnt, sendtype, recvbuf, recvcount, recvtype, root, comm )
void * sendbuf;
int sendcnt;
MPI_Datatype sendtype;
void * recvbuf;
int recvcount;
MPI_Datatype recvtype;
int root;
MPI_Comm comm;
{
    int np, nb; 
    PMPI_Comm_size( comm, &np);
    nb = E_count2byte( sendtype, sendcnt);
    E_GET_COLL_AVG_TIME(n_gather,gather)
}
double E_MPI_Gatherv( sendbuf, sendcnt, sendtype, recvbuf, recvcnts, displs, recvtype, root, comm )
void * sendbuf;
int sendcnt;
MPI_Datatype sendtype;
void * recvbuf;
int * recvcnts;
int * displs;
MPI_Datatype recvtype;
int root;
MPI_Comm comm;
{
    int np, nb; 
    PMPI_Comm_size( comm, &np);
    nb = E_count2byte( sendtype, sendcnt);
    E_GET_COLL_AVG_TIME(n_gatherv,gatherv)
}
double E_MPI_Reduce_scatter( sendbuf, recvbuf, recvcnts, datatype, op, comm )
void * sendbuf;
void * recvbuf;
int * recvcnts;
MPI_Datatype datatype;
MPI_Op op;
MPI_Comm comm;
{
    return 0;
}
double E_MPI_Reduce( sendbuf, recvbuf, count, datatype, op, root, comm )
void * sendbuf;
void * recvbuf;
int count;
MPI_Datatype datatype;
MPI_Op op;
int root;
MPI_Comm comm;
{
    int np, nb; 
    PMPI_Comm_size( comm, &np);
    nb = E_count2byte( datatype, count );
    E_GET_COLL_AVG_TIME(n_reduce,reduce)
}
double E_MPI_Scan( sendbuf, recvbuf, count, datatype, op, comm )
void * sendbuf;
void * recvbuf;
int count;
MPI_Datatype datatype;
MPI_Op op;
MPI_Comm comm;
{
	return 0;
}
double E_MPI_Scatter( sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm )
void * sendbuf;
int sendcnt;
MPI_Datatype sendtype;
void * recvbuf;
int recvcnt;
MPI_Datatype recvtype;
int root;
MPI_Comm comm;
{
    int np, nb; 
    PMPI_Comm_size( comm, &np);
    nb = E_count2byte( sendtype, sendcnt );
    E_GET_COLL_AVG_TIME(n_scatter,scatter);
}
double E_MPI_Scatterv( sendbuf, sendcnts, displs, sendtype, recvbuf, recvcnt, recvtype, root, comm )
void * sendbuf;
int * sendcnts;
int * displs;
MPI_Datatype sendtype;
void * recvbuf;
int recvcnt;
MPI_Datatype recvtype;
int root;
MPI_Comm comm;
{
    int np, nb; 
    PMPI_Comm_size( comm, &np);
    nb = E_count2byte( sendtype, *sendcnts );
    E_GET_COLL_AVG_TIME(n_scatterv,scatterv);
}
/******************************************************************
*                                                                 *
*               MPI Communicator Related Functions                *
*                                                                 *
******************************************************************/
double E_MPI_Comm_compare( comm1, comm2, result )
MPI_Comm comm1;
MPI_Comm comm2;
int * result;
{
	return 0;
}
double E_MPI_Comm_create( comm, group, comm_out )
MPI_Comm comm;
MPI_Group group;
MPI_Comm * comm_out;
{
	return 0;
}
double E_MPI_Comm_dup( comm, comm_out )
MPI_Comm comm;
MPI_Comm * comm_out;
{
	return 0;
}
double E_MPI_Comm_free( comm )
MPI_Comm * comm;
{
	return 0;
}
double E_MPI_Comm_group( comm, group )
MPI_Comm comm;
MPI_Group * group;
{
	return 0;
}
double E_MPI_Comm_rank( comm, rank )
MPI_Comm comm;
int * rank;
{
	return 0;
}
double E_MPI_Comm_remote_group( comm, group )
MPI_Comm comm;
MPI_Group * group;
{
	return 0;
}
double E_MPI_Comm_remote_size( comm, size )
MPI_Comm comm;
int * size;
{
	return 0;
}
double E_MPI_Comm_size( comm, size )
MPI_Comm comm;
int * size;
{
	return 0;
}
double E_MPI_Comm_split( comm, color, key, comm_out )
MPI_Comm comm;
int color;
int key;
MPI_Comm * comm_out;
{
	return 0;
}
double E_MPI_Comm_test_inter( comm, flag )
MPI_Comm comm;
int * flag;
{
	return 0;
}
double E_MPI_Group_compare( group1, group2, result )
MPI_Group group1;
MPI_Group group2;
int * result;
{
	return 0;
}
double E_MPI_Group_difference( group1, group2, group_out )
MPI_Group group1;
MPI_Group group2;
MPI_Group * group_out;
{
	return 0;
}
double E_MPI_Group_excl( group, n, ranks, newgroup )
MPI_Group group;
int n;
int * ranks;
MPI_Group * newgroup;
{
	return 0;
}
double E_MPI_Group_free( group )
MPI_Group * group;
{
	return 0;
}
double E_MPI_Group_incl( group, n, ranks, group_out )
MPI_Group group;
int n;
int * ranks;
MPI_Group * group_out;
{
	return 0;
}
double E_MPI_Group_intersection( group1, group2, group_out )
MPI_Group group1;
MPI_Group group2;
MPI_Group * group_out;
{
	return 0;
}
double E_MPI_Group_rank( group, rank )
MPI_Group group;
int * rank;
{
	return 0;
}
double E_MPI_Group_range_excl( group, n, ranges, newgroup )
MPI_Group group;
int n;
int ranges[][3];
MPI_Group * newgroup;
{
	return 0;
}
double E_MPI_Group_range_incl( group, n, ranges, newgroup )
MPI_Group group;
int n;
int ranges[][3];
MPI_Group * newgroup;
{
	return 0;
}
double E_MPI_Group_size( group, size )
MPI_Group group;
int * size;
{
	return 0;
}
double E_MPI_Group_translate_ranks( group_a, n, ranks_a, group_b, ranks_b )
MPI_Group group_a;
int n;
int * ranks_a;
MPI_Group group_b;
int * ranks_b;
{
	return 0;
}
double E_MPI_Group_union( group1, group2, group_out )
MPI_Group group1;
MPI_Group group2;
MPI_Group * group_out;
{
	return 0;
}
double E_MPI_Intercomm_create( local_comm, local_leader, peer_comm, remote_leader, tag, comm_out )
MPI_Comm local_comm;
int local_leader;
MPI_Comm peer_comm;
int remote_leader;
int tag;
MPI_Comm * comm_out;
{
	return 0;
}
double E_MPI_Intercomm_merge( comm, high, comm_out )
MPI_Comm comm;
int high;
MPI_Comm * comm_out;
{
	return 0;
}
double E_MPI_Keyval_create( copy_fn, delete_fn, keyval, extra_state )
MPI_Copy_function * copy_fn;
MPI_Delete_function * delete_fn;
int * keyval;
void * extra_state;
{
	return 0;
}
double E_MPI_Keyval_free( keyval )
int * keyval;
{
	return 0;
}
double E_MPI_Cart_coords( comm, rank, maxdims, coords )
MPI_Comm comm;
int rank;
int maxdims;
int * coords;
{
	return 0;
}
double E_MPI_Cart_create( comm_old, ndims, dims, periods, reorder, comm_cart )
MPI_Comm comm_old;
int ndims;
int * dims;
int * periods;
int reorder;
MPI_Comm * comm_cart;
{
	return 0;
}
double E_MPI_Cart_get( comm, maxdims, dims, periods, coords )
MPI_Comm comm;
int maxdims;
int * dims;
int * periods;
int * coords;
{
	return 0;
}
double E_MPI_Cart_map( comm_old, ndims, dims, periods, newrank )
MPI_Comm comm_old;
int ndims;
int * dims;
int * periods;
int * newrank;
{
	return 0;
}
double E_MPI_Cart_rank( comm, coords, rank )
MPI_Comm comm;
int * coords;
int * rank;
{
	return 0;
}
double E_MPI_Cart_shift( comm, direction, displ, source, dest )
MPI_Comm comm;
int direction;
int displ;
int * source;
int * dest;
{
	return 0;
}
double E_MPI_Cart_sub( comm, remain_dims, comm_new )
MPI_Comm comm;
int * remain_dims;
MPI_Comm * comm_new;
{
	return 0;
}
double E_MPI_Cartdim_get( comm, ndims )
MPI_Comm comm;
int * ndims;
{
	return 0;
}
double E_MPI_Dims_create( nnodes, ndims, dims )
int nnodes;
int ndims;
int * dims;
{
	return 0;
}
double E_MPI_Graph_create( comm_old, nnodes, index, edges, reorder, comm_graph )
MPI_Comm comm_old;
int nnodes;
int * index;
int * edges;
int reorder;
MPI_Comm * comm_graph;
{
	return 0;
}
double E_MPI_Graph_get( comm, maxindex, maxedges, index, edges )
MPI_Comm comm;
int maxindex;
int maxedges;
int * index;
int * edges;
{
	return 0;
}
double E_MPI_Graph_map( comm_old, nnodes, index, edges, newrank )
MPI_Comm comm_old;
int nnodes;
int * index;
int * edges;
int * newrank;
{
	return 0;
}
double E_MPI_Graph_neighbors( comm, rank, maxneighbors, neighbors )
MPI_Comm comm;
int rank;
int  maxneighbors;
int * neighbors;
{
	return 0;
}
double E_MPI_Graph_neighbors_count( comm, rank, nneighbors )
MPI_Comm comm;
int rank;
int * nneighbors;
{
	return 0;
}
double E_MPI_Graphdims_get( comm, nnodes, nedges )
MPI_Comm comm;
int * nnodes;
int * nedges;
{
	return 0;
}
double E_MPI_Topo_test( comm, top_type )
MPI_Comm comm;
int * top_type;
{
	return 0;
}
/******************************************************************
*                                                                 *
*                      MPI Other Functions                        *
*                                                                 *
******************************************************************/
/* LAM MPI defines MPI_Abort as a macro! We check for this and if 
 */
#if (defined(MPI_Abort) && defined(_ULM_MPI_H_))
int E_MPI_Abort( MPI_Comm comm, int errorcode, char * file, int line)
#else
int E_MPI_Abort( comm, errorcode )
MPI_Comm comm;
int errorcode;
#endif /* MPI_Abort & LAM MPI [LAM MPI] */
{
	return 0;
}
double E_MPI_Error_class( errorcode, errorclass )
int errorcode;
int * errorclass;
{
	return 0;
}
double E_MPI_Errhandler_create( function, errhandler )
MPI_Handler_function * function;
MPI_Errhandler * errhandler;
{
	return 0;
}
double E_MPI_Errhandler_free( errhandler )
MPI_Errhandler * errhandler;
{
	return 0;
}
double E_MPI_Errhandler_get( comm, errhandler )
MPI_Comm comm;
MPI_Errhandler * errhandler;
{
	return 0;
}
double E_MPI_Error_string( errorcode, string, resultlen )
int errorcode;
char * string;
int * resultlen;
{
	return 0;
}
double E_MPI_Errhandler_set( comm, errhandler )
MPI_Comm comm;
MPI_Errhandler errhandler;
{
	return 0;
}
double E_MPI_Get_processor_name( name, resultlen )
char * name;
int * resultlen;
{
	return 0;
}
double  E_MPI_Wtick(  )
{
	return 0;
}
double  E_MPI_Wtime(  )
{
	return 0;
}
double E_MPI_Address( location, address )
void * location;
MPI_Aint * address;
{
	return 0;
}
double E_MPI_Op_create( function, commute, op )
MPI_User_function * function;
int commute;
MPI_Op * op;
{
	return 0;
}
double E_MPI_Op_free( op )
MPI_Op * op;
{
	return 0;
}
double E_MPI_Attr_delete( comm, keyval )
MPI_Comm comm;
int keyval;
{
	return 0;
}
double E_MPI_Attr_get( comm, keyval, attr_value, flag )
MPI_Comm comm;
int keyval;
void * attr_value;
int * flag;
{
	return 0;
}
double E_MPI_Attr_put( comm, keyval, attr_value )
MPI_Comm comm;
int keyval;
void * attr_value;
{
	return 0;
}
double E_MPI_Buffer_attach( buffer, size )
void * buffer;
int size;
{
	return 0;
}
double E_MPI_Buffer_detach( buffer, size )
void * buffer;
int * size;
{
	return 0;
}
double E_MPI_Get_elements( status, datatype, elements )
MPI_Status * status;
MPI_Datatype datatype;
int * elements;
{
	return 0;
}
double E_MPI_Get_count( status, datatype, count )
MPI_Status * status;
MPI_Datatype datatype;
int * count;
{
	return 0;
}
double E_MPI_Type_commit( datatype )
MPI_Datatype * datatype;
{
	return 0;
}
double E_MPI_Type_contiguous( count, old_type, newtype )
int count;
MPI_Datatype old_type;
MPI_Datatype * newtype;
{
	return 0;
}
double E_MPI_Type_extent( datatype, extent )
MPI_Datatype datatype;
MPI_Aint * extent;
{
	return 0;
}
double E_MPI_Type_free( datatype )
MPI_Datatype * datatype;
{
	return 0;
}
double E_MPI_Type_hindexed( count, blocklens, indices, old_type, newtype )
int count;
int * blocklens;
MPI_Aint * indices;
MPI_Datatype old_type;
MPI_Datatype * newtype;
{
	return 0;
}
double E_MPI_Type_hvector( count, blocklen, stride, old_type, newtype )
int count;
int blocklen;
MPI_Aint stride;
MPI_Datatype old_type;
MPI_Datatype * newtype;
{
	return 0;
}
double E_MPI_Type_indexed( count, blocklens, indices, old_type, newtype )
int count;
int * blocklens;
int * indices;
MPI_Datatype old_type;
MPI_Datatype * newtype;
{
	return 0;
	return 0;
}
double E_MPI_Type_lb( datatype, displacement )
MPI_Datatype datatype;
MPI_Aint * displacement;
{
	return 0;
	return 0;
}
double E_MPI_Type_size( datatype, size )
MPI_Datatype datatype;
int * size;
{
	return 0;
}
double E_MPI_Type_struct( count, blocklens, indices, old_types, newtype )
int count;
int * blocklens;
MPI_Aint * indices;
MPI_Datatype * old_types;
MPI_Datatype * newtype;
{
	return 0;
}
double E_MPI_Type_ub( datatype, displacement )
MPI_Datatype datatype;
MPI_Aint * displacement;
{
	return 0;
}
double E_MPI_Type_vector( count, blocklen, stride, old_type, newtype )
int count;
int blocklen;
int stride;
MPI_Datatype old_type;
MPI_Datatype * newtype;
{
	return 0;
}
double E_MPI_Unpack( inbuf, insize, position, outbuf, outcount, type, comm )
void * inbuf;
int insize;
int * position;
void * outbuf;
int outcount;
MPI_Datatype type;
MPI_Comm comm;
{
	return 0;
}
double E_MPI_Pack( inbuf, incount, type, outbuf, outcount, position, comm )
void * inbuf;
int incount;
MPI_Datatype type;
void * outbuf;
int outcount;
int * position;
MPI_Comm comm;
{
	return 0;
}
double E_MPI_Pack_size( incount, datatype, comm, size )
int incount;
MPI_Datatype datatype;
MPI_Comm comm;
int * size;
{
	return 0;
}
/********************************************************
*                  User defined functions               *
********************************************************/
void E_MPI_Profile_on()
{
}
void E_MPI_Profile_off()
{
}
