/********************************************************************
 *                                                                   *  
 * This library can be used to acquire communication and computation *
 * profile for a given MPI program. Written by Jidong Zhai.          *
 * July.12.2012. Contact: zhaijidong@gmail.com                       *
 *                                                                   *
 ********************************************************************/

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "timer.h"
#include "estimator.h"
#include "reporter.h"

/******************************************************************
 *                                                                 *
 *                   MPI Functions for Management                  *
 *                                                                 *
 ******************************************************************/

int  MPI_Init( argc, argv )
    int * argc;
    char *** argv;
{
    int  returnVal;
    int  proc_id;
    returnVal = PMPI_Init( argc, argv );

#ifdef PERF_ASSERT
    PMPI_Comm_rank( MPI_COMM_WORLD, &proc_id );
    E_MPI_Init(argc, argv);
    char filename[128];
    sprintf(filename,"%s%d","mpipa.report.",proc_id);
    R_init(filename);
#endif
    return returnVal;
}

#ifdef PERF_MPI_THREADED
int  MPI_Init_thread (argc, argv, required, provided ) //1
    int * argc;
    char *** argv;
    int required;
    int *provided;
{
    int  returnVal;
    returnVal = PMPI_Init_thread( argc, argv, required, provided );
    return returnVal;
}
#endif /* PERF_MPI_THREADED */


int  MPI_Finalize(  )
{
    int  returnVal;

#ifdef PERF_ASSERT
    E_MPI_Finalize();
    int numproc = 0;
    PMPI_Comm_size(MPI_COMM_WORLD, &numproc);
    R_report(R_Level,numproc);
#endif

    returnVal = PMPI_Finalize();

    return returnVal;
}


/******************************************************************
 *                                                                 *
 *          MPI Point-to-Point(P2P) Communication Related          *
 *                                                                 *
 ******************************************************************/

int  MPI_Bsend( buf, count, datatype, dest, tag, comm ) //3
    void * buf;
    int count;
    MPI_Datatype datatype;
    int dest;
    int tag;
    MPI_Comm comm;
{
    int returnVal;
    returnVal = PMPI_Bsend( buf, count, datatype, dest, tag, comm );
    return returnVal;
}

int  MPI_Bsend_init( buf, count, datatype, dest, tag, comm, request ) //4
    void * buf;
    int count;
    MPI_Datatype datatype;
    int dest;
    int tag;
    MPI_Comm comm;
    MPI_Request * request;
{
    int  returnVal;
    returnVal = PMPI_Bsend_init( buf, count, datatype, dest, tag, comm, request );
    return returnVal;
}

int  MPI_Recv_init( buf, count, datatype, source, tag, comm, request ) //5
    void * buf;
    int count;
    MPI_Datatype datatype;
    int source;
    int tag;
    MPI_Comm comm;
    MPI_Request * request;
{
    int  returnVal;
    returnVal = PMPI_Recv_init( buf, count, datatype, source, tag, comm, request );
    return returnVal;
}

int  MPI_Send_init( buf, count, datatype, dest, tag, comm, request ) //6
    void * buf;
    int count;
    MPI_Datatype datatype;
    int dest;
    int tag;
    MPI_Comm comm;
    MPI_Request * request;
{
    int  returnVal;
    returnVal = PMPI_Send_init( buf, count, datatype, dest, tag, comm, request );
    return returnVal;
}

int  MPI_Ibsend( buf, count, datatype, dest, tag, comm, request ) //7
    void * buf;
    int count;
    MPI_Datatype datatype;
    int dest;
    int tag;
    MPI_Comm comm;
    MPI_Request * request;
{
    int  returnVal;
    returnVal = PMPI_Ibsend( buf, count, datatype, dest, tag, comm, request );
    return returnVal;
}

int  MPI_Irecv( buf, count, datatype, source, tag, comm, request ) //8
    void * buf;
    int count;
    MPI_Datatype datatype;
    int source;
    int tag;
    MPI_Comm comm;
    MPI_Request * request;
{
    int  returnVal;
#ifdef PERF_ASSERT
    _timer_start(PATN);
#endif 
    returnVal = PMPI_Irecv( buf, count, datatype, source, tag, comm, request );
#ifdef PERF_ASSERT
    _timer_stop(PATN);
    double r = _timer_read(PATN);
    _timer_clear(PATN);
    double e = E_MPI_Irecv( buf, count, datatype, source, tag, comm, request );
    int pid;
    PMPI_Comm_rank(MPI_COMM_WORLD, &pid);
    char para[100];
    sprintf(para, "count:%d, source:%d",count,source);
    R_log(R_Level,0,r,e,pid,8,para);
#endif
    return returnVal;
}

int  MPI_Irsend( buf, count, datatype, dest, tag, comm, request ) //9
    void * buf;
    int count;
    MPI_Datatype datatype;
    int dest;
    int tag;
    MPI_Comm comm;
    MPI_Request * request;
{
    int  returnVal;
    returnVal = PMPI_Irsend( buf, count, datatype, dest, tag, comm, request );
    return returnVal;
}

int  MPI_Isend( buf, count, datatype, dest, tag, comm, request ) //10
    void * buf;
    int count;
    MPI_Datatype datatype;
    int dest;
    int tag;
    MPI_Comm comm;
    MPI_Request * request;
{
    int  returnVal;
#ifdef PERF_ASSERT
    _timer_start(PATN);
#endif 
    returnVal = PMPI_Isend( buf, count, datatype, dest, tag, comm, request );
#ifdef PERF_ASSERT
    _timer_stop(PATN);
    double r = _timer_read(PATN);
    _timer_clear(PATN);
    double e = E_MPI_Isend( buf, count, datatype, dest, tag, comm, request );
    int pid;
    PMPI_Comm_rank(MPI_COMM_WORLD, &pid);
    char para[100];
    sprintf(para, "count:%d, dest:%d",count,dest);
    R_log(R_Level,0,r,e,pid,10,para);
#endif
    return returnVal;
}

int  MPI_Issend( buf, count, datatype, dest, tag, comm, request ) //11
    void * buf;
    int count;
    MPI_Datatype datatype;
    int dest;
    int tag;
    MPI_Comm comm;
    MPI_Request * request;
{
    int  returnVal;
    returnVal = PMPI_Issend( buf, count, datatype, dest, tag, comm, request );
    return returnVal;
}

int  MPI_Recv( buf, count, datatype, source, tag, comm, status ) //12
    void * buf;
    int count;
    MPI_Datatype datatype;
    int source;
    int tag;
    MPI_Comm comm;
    MPI_Status * status;
{
    int  returnVal;
#ifdef PERF_ASSERT
    _timer_start(PATN);
#endif 
    returnVal = PMPI_Recv( buf, count, datatype, source, tag, comm, status );
#ifdef PERF_ASSERT
    _timer_stop(PATN);
    double r = _timer_read(PATN);
    _timer_clear(PATN);
    double e = E_MPI_Recv( buf, count, datatype, source, tag, comm, status );
    int pid;
    PMPI_Comm_rank(MPI_COMM_WORLD, &pid);
    char para[100];
    sprintf(para, "count:%d, source:%d",count,source);
    R_log(R_Level,0,r,e,pid,12,para);
#endif
    return returnVal;
}

int  MPI_Rsend( buf, count, datatype, dest, tag, comm ) //13
    void * buf;
    int count;
    MPI_Datatype datatype;
    int dest;
    int tag;
    MPI_Comm comm;
{
    int  returnVal;
    returnVal = PMPI_Rsend( buf, count, datatype, dest, tag, comm );
    return returnVal;
}

int  MPI_Rsend_init( buf, count, datatype, dest, tag, comm, request ) //14
    void * buf;
    int count;
    MPI_Datatype datatype;
    int dest;
    int tag;
    MPI_Comm comm;
    MPI_Request * request;
{
    int  returnVal;
    returnVal = PMPI_Rsend_init( buf, count, datatype, dest, tag, comm, request );
    return returnVal;
}


int  MPI_Send( buf, count, datatype, dest, tag, comm ) //15
    void * buf;
    int count;
    MPI_Datatype datatype;
    int dest;
    int tag;
    MPI_Comm comm;
{
    int  returnVal;
#ifdef PERF_ASSERT
    _timer_start(PATN);
#endif 
    returnVal = PMPI_Send( buf, count, datatype, dest, tag, comm );
#ifdef PERF_ASSERT
    _timer_stop(PATN);
    double r = _timer_read(PATN);
    _timer_clear(PATN);
    double e = E_MPI_Send( buf, count, datatype, dest, tag, comm );
    int pid;
    PMPI_Comm_rank(MPI_COMM_WORLD, &pid);
    char para[100];
    sprintf(para, "count:%d, dest:%d",count,dest);
    R_log(R_Level,0,r,e,pid,15,para);
#endif
    return returnVal;
}

int  MPI_Sendrecv( sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status ) //16
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
    int  returnVal;
#ifdef PERF_ASSERT
    _timer_start(PATN);
#endif 
    returnVal = PMPI_Sendrecv( sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status );
#ifdef PERF_ASSERT
    _timer_stop(PATN);
    double r = _timer_read(PATN);
    _timer_clear(PATN);
    double e = E_MPI_Sendrecv( sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status );
    int pid;
    PMPI_Comm_rank(MPI_COMM_WORLD, &pid);
    char para[100];
    sprintf(para, "sendcount:%d, dest:%d, recvcount:%d, source:%d",sendcount,dest,recvcount,source);
    R_log(R_Level,0,r,e,pid,16,para);
#endif
    return returnVal;
}

int  MPI_Sendrecv_replace( buf, count, datatype, dest, sendtag, source, recvtag, comm, status ) //17
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
    int  returnVal;
    returnVal = PMPI_Sendrecv_replace( buf, count, datatype, dest, sendtag, source, recvtag, comm, status );
    return returnVal;
}

int  MPI_Ssend( buf, count, datatype, dest, tag, comm ) //18
    void * buf;
    int count;
    MPI_Datatype datatype;
    int dest;
    int tag;
    MPI_Comm comm;
{
    int  returnVal;
    returnVal = PMPI_Ssend( buf, count, datatype, dest, tag, comm );
    return returnVal;
}

int  MPI_Ssend_init( buf, count, datatype, dest, tag, comm, request ) //19
    void * buf;
    int count;
    MPI_Datatype datatype;
    int dest;
    int tag;
    MPI_Comm comm;
    MPI_Request * request;
{
    int  returnVal;
    returnVal = PMPI_Ssend_init( buf, count, datatype, dest, tag, comm, request );
    return returnVal;
}


int   MPI_Test( request, flag, status ) //20
    MPI_Request * request;
    int * flag;
    MPI_Status * status;
{
    int   returnVal;
    returnVal = PMPI_Test( request, flag, status );
    return returnVal;
}

int  MPI_Testall( count, array_of_requests, flag, array_of_statuses ) //21
    int count;
    MPI_Request * array_of_requests;
    int * flag;
    MPI_Status * array_of_statuses;
{
    int returnVal;
    returnVal = PMPI_Testall( count, array_of_requests, flag, array_of_statuses );
    return returnVal;
}

int  MPI_Testany( count, array_of_requests, index, flag, status ) //22
    int count;
    MPI_Request * array_of_requests;
    int * index;
    int * flag;
    MPI_Status * status;
{
    int  returnVal;
    returnVal = PMPI_Testany( count, array_of_requests, index, flag, status );
    return returnVal;
}

int  MPI_Test_cancelled( status, flag ) //23
    MPI_Status * status;
    int * flag;
{
    int  returnVal;
    returnVal = PMPI_Test_cancelled( status, flag );
    return returnVal;
}

int  MPI_Testsome( incount, array_of_requests, outcount, array_of_indices, array_of_statuses ) //24
    int incount;
    MPI_Request * array_of_requests;
    int * outcount;
    int * array_of_indices;
    MPI_Status * array_of_statuses;
{
    int  returnVal;
    returnVal = PMPI_Testsome( incount, array_of_requests, outcount, array_of_indices, array_of_statuses );
    return returnVal;
}

int   MPI_Wait( request, status ) //25
    MPI_Request * request;
    MPI_Status * status;
{
    int   returnVal;
    //MPI_Status local_status;
    //MPI_Request* saverequest = request;

#ifdef PERF_ASSERT
    double e = E_MPI_Wait( request, status );
    _timer_start(PATN);
#endif 
    returnVal = PMPI_Wait( request, status );
#ifdef PERF_ASSERT
    _timer_stop(PATN);
    double r = _timer_read(PATN);
    _timer_clear(PATN);
    int pid;
    PMPI_Comm_rank(MPI_COMM_WORLD, &pid);
    char para[100];
    sprintf(para, "req:%p", request);
    R_log(R_Level,0,r,e,pid,25,para);
#endif
    return returnVal;
}

int  MPI_Waitall( count, array_of_requests, array_of_statuses ) //26
    int count;
    MPI_Request * array_of_requests;
    MPI_Status * array_of_statuses;
{
    int  returnVal;
#ifdef PERF_ASSERT
    double e = E_MPI_Waitall( count, array_of_requests, array_of_statuses );
    _timer_start(PATN);
#endif 
    returnVal = PMPI_Waitall( count, array_of_requests, array_of_statuses );
#ifdef PERF_ASSERT
    _timer_stop(PATN);
    double r = _timer_read(PATN);
    _timer_clear(PATN);
    int pid;
    PMPI_Comm_rank(MPI_COMM_WORLD, &pid);
    char para[100];
    sprintf(para, "count:%d", count);
    R_log(R_Level,0,r,e,pid,26,para);
#endif
    return returnVal;
}

int  MPI_Waitany( count, array_of_requests, index, status ) //27
    int count;
    MPI_Request * array_of_requests;
    int * index;
    MPI_Status * status;
{
    int  returnVal;
    returnVal = PMPI_Waitany( count, array_of_requests, index, status );
    return returnVal;
}

int  MPI_Waitsome( incount, array_of_requests, outcount, array_of_indices, array_of_statuses ) //28
    int incount;
    MPI_Request * array_of_requests;
    int * outcount;
    int * array_of_indices;
    MPI_Status * array_of_statuses;
{
    int  returnVal;
    returnVal = PMPI_Waitsome( incount, array_of_requests, outcount, array_of_indices, array_of_statuses );
    return returnVal;
}

int  MPI_Cancel( request ) //29
    MPI_Request * request;
{
    int  returnVal;
    returnVal = PMPI_Cancel( request );
    return returnVal;
}

int  MPI_Request_free( request ) //30
    MPI_Request * request;
{
    int  returnVal;
    returnVal = PMPI_Request_free( request );
    return returnVal;
}

int  MPI_Start( request ) //31
    MPI_Request * request;
{
    int  returnVal;
    returnVal = PMPI_Start( request );
    return returnVal;
}

int  MPI_Startall( count, array_of_requests ) //32
    int count;
    MPI_Request * array_of_requests;
{
    int  returnVal;
    returnVal = PMPI_Startall( count, array_of_requests );
    return returnVal;
}

int  MPI_Iprobe( source, tag, comm, flag, status ) //33
    int source;
    int tag;
    MPI_Comm comm;
    int * flag;
    MPI_Status * status;
{
    int  returnVal;
    returnVal = PMPI_Iprobe( source, tag, comm, flag, status );
    return returnVal;
}

int  MPI_Probe( source, tag, comm, status ) //34
    int source;
    int tag;
    MPI_Comm comm;
    MPI_Status * status;
{
    int  returnVal;
    returnVal = PMPI_Probe( source, tag, comm, status );
    return returnVal;
}

/******************************************************************
 *                                                                 *
 *               MPI Collective Communication Related              *
 *                                                                 *
 ******************************************************************/

int   MPI_Allgather( sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm ) //35
    void * sendbuf;
    int sendcount;
    MPI_Datatype sendtype;
    void * recvbuf;
    int recvcount;
    MPI_Datatype recvtype;
    MPI_Comm comm;
{
    int   returnVal;

#ifdef PERF_ASSERT
    _timer_start(PATN);
#endif 

    returnVal = PMPI_Allgather( sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm );

#ifdef PERF_ASSERT
    _timer_stop(PATN);
#endif 

#ifdef PERF_ASSERT
    double r = _timer_read(PATN);
    _timer_clear(PATN);
    double e = E_MPI_Allgather( sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm );
    int pid;
    PMPI_Comm_rank(MPI_COMM_WORLD, &pid);
    char para[100];
    sprintf(para, "sendcount:%d, recvcount:%d",sendcount,recvcount);
    R_log(R_Level,0,r,e,pid,35,para);
#endif

    return returnVal;
}

int   MPI_Allgatherv( sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm ) //36
    void * sendbuf;
    int sendcount;
    MPI_Datatype sendtype;
    void * recvbuf;
    int * recvcounts;
    int * displs;
    MPI_Datatype recvtype;
    MPI_Comm comm;
{
    int   returnVal;

#ifdef PERF_ASSERT
    _timer_start(PATN);
#endif 

    returnVal = PMPI_Allgatherv( sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm );

#ifdef PERF_ASSERT
    _timer_stop(PATN);
#endif 

#ifdef PERF_ASSERT
    double r = _timer_read(PATN);
    _timer_clear(PATN);
    double e = E_MPI_Allgatherv( sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm );
    int pid;
    PMPI_Comm_rank(MPI_COMM_WORLD, &pid);
    char para[100];
    sprintf(para, "sendcount:%d, recvcounts:%d",sendcount,*recvcounts);
    R_log(R_Level,0,r,e,pid,36,para);
#endif

    return returnVal;
}

int   MPI_Allreduce( sendbuf, recvbuf, count, datatype, op, comm ) //37
    void * sendbuf;
    void * recvbuf;
    int count;
    MPI_Datatype datatype;
    MPI_Op op;
    MPI_Comm comm;
{
    int   returnVal;

#ifdef PERF_ASSERT
    _timer_start(PATN);
#endif 

    returnVal = PMPI_Allreduce( sendbuf, recvbuf, count, datatype, op, comm );

#ifdef PERF_ASSERT
    _timer_stop(PATN);
#endif 

#ifdef PERF_ASSERT
    double r = _timer_read(PATN);
    _timer_clear(PATN);
    double e = E_MPI_Allreduce( sendbuf, recvbuf, count, datatype, op, comm );
    int pid;
    PMPI_Comm_rank(MPI_COMM_WORLD, &pid);
    char para[100];
    sprintf(para, "count:%d",count);
    R_log(R_Level,0,r,e,pid,37,para);
#endif

    return returnVal;
}

int  MPI_Alltoall( sendbuf, sendcount, sendtype, recvbuf, recvcnt, recvtype, comm ) //38
    void * sendbuf;
    int sendcount;
    MPI_Datatype sendtype;
    void * recvbuf;
    int recvcnt;
    MPI_Datatype recvtype;
    MPI_Comm comm;
{
    int  returnVal;

#ifdef PERF_ASSERT
    _timer_start(PATN);
#endif 

    returnVal = PMPI_Alltoall( sendbuf, sendcount, sendtype, recvbuf, recvcnt, recvtype, comm );

#ifdef PERF_ASSERT
    _timer_stop(PATN);
#endif 

#ifdef PERF_ASSERT
    double r = _timer_read(PATN);
    _timer_clear(PATN);
    double e = E_MPI_Alltoall( sendbuf, sendcount, sendtype, recvbuf, recvcnt, recvtype, comm );
    int pid;
    PMPI_Comm_rank(MPI_COMM_WORLD, &pid);
    char para[100];
    sprintf(para, "sendcount:%d, recvcnt:%d",sendcount, recvcnt);
    R_log(R_Level,0,r,e,pid,38,para);
#endif

    return returnVal;
}

int   MPI_Alltoallv( sendbuf, sendcnts, sdispls, sendtype, recvbuf, recvcnts, rdispls, recvtype, comm ) //39
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
    int   returnVal;

#ifdef PERF_ASSERT
    _timer_start(PATN);
#endif 

    returnVal = PMPI_Alltoallv( sendbuf, sendcnts, sdispls, sendtype, recvbuf, recvcnts, rdispls, recvtype, comm );

#ifdef PERF_ASSERT
    _timer_stop(PATN);
#endif 

#ifdef PERF_ASSERT
    double r = _timer_read(PATN);
    _timer_clear(PATN);
    double e = E_MPI_Alltoallv( sendbuf, sendcnts, sdispls, sendtype, recvbuf, recvcnts, rdispls, recvtype, comm );
    int pid;
    PMPI_Comm_rank(MPI_COMM_WORLD, &pid);
    char para[100];
    sprintf(para, "sendcnts:%d, recvcnts:%d",*sendcnts, *recvcnts);
    R_log(R_Level,0,r,e,pid,39,para);
#endif

    return returnVal;
}

int   MPI_Barrier( comm ) //40
    MPI_Comm comm;
{
    int   returnVal;

#ifdef PERF_ASSERT
    _timer_start(PATN);
#endif 

    returnVal = PMPI_Barrier( comm );

#ifdef PERF_ASSERT
    _timer_stop(PATN);
#endif 

#ifdef PERF_ASSERT
    double r = _timer_read(PATN);
    _timer_clear(PATN);
    double e = E_MPI_Barrier( comm );
    int pid;
    PMPI_Comm_rank(MPI_COMM_WORLD, &pid);
    char para[100];
    int size;
    PMPI_Comm_size(comm,&size);
    sprintf(para, "CommSize:%d",size);
    R_log(R_Level,0,r,e,pid,40,para);
#endif

    return returnVal;
}

int   MPI_Bcast( buffer, count, datatype, root, comm ) //41
    void * buffer;
    int count;
    MPI_Datatype datatype;
    int root;
    MPI_Comm comm;
{
    int   returnVal;

#ifdef PERF_ASSERT
    _timer_start(PATN);
#endif 
    returnVal = PMPI_Bcast( buffer, count, datatype, root, comm );

#ifdef PERF_ASSERT
    _timer_stop(PATN);
    double r = _timer_read(PATN);
    _timer_clear(PATN);
    double e = E_MPI_Bcast( buffer, count, datatype, root, comm );
    int pid;
    PMPI_Comm_rank(MPI_COMM_WORLD, &pid);
    char para[100];
    sprintf(para, "count:%d, root:%d",count,root);
    R_log(R_Level,0,r,e,pid,41,para);
#endif

    return returnVal;
}

int   MPI_Gather( sendbuf, sendcnt, sendtype, recvbuf, recvcount, recvtype, root, comm ) //42
    void * sendbuf;
    int sendcnt;
    MPI_Datatype sendtype;
    void * recvbuf;
    int recvcount;
    MPI_Datatype recvtype;
    int root;
    MPI_Comm comm;
{
    int   returnVal;

#ifdef PERF_ASSERT
    _timer_start(PATN);
#endif 

    returnVal = PMPI_Gather( sendbuf, sendcnt, sendtype, recvbuf, recvcount, recvtype, root, comm );

#ifdef PERF_ASSERT
    _timer_stop(PATN);
    double r = _timer_read(PATN);
    _timer_clear(PATN);
    double e = E_MPI_Gather( sendbuf, sendcnt, sendtype, recvbuf, recvcount, recvtype, root, comm );
    int pid;
    PMPI_Comm_rank(MPI_COMM_WORLD, &pid);
    char para[100];
    sprintf(para, "sendcnt:%d, recvcount:%d, root:%d",sendcnt,recvcount,root);
    R_log(R_Level,0,r,e,pid,42,para);
#endif

    return returnVal;
}

int   MPI_Gatherv( sendbuf, sendcnt, sendtype, recvbuf, recvcnts, displs, recvtype, root, comm ) //43
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
    int   returnVal;

#ifdef PERF_ASSERT
    _timer_start(PATN);
#endif 

    returnVal = PMPI_Gatherv( sendbuf, sendcnt, sendtype, recvbuf, recvcnts, displs, recvtype, root, comm );

#ifdef PERF_ASSERT
    _timer_stop(PATN);
    double r = _timer_read(PATN);
    _timer_clear(PATN);
    double e = E_MPI_Gatherv( sendbuf, sendcnt, sendtype, recvbuf, recvcnts, displs, recvtype, root, comm );
    int pid;
    PMPI_Comm_rank(MPI_COMM_WORLD, &pid);
    char para[100];
    sprintf(para, "sendcnt:%d, recvcnts:%d, root:%d",sendcnt,*recvcnts,root);
    R_log(R_Level,0,r,e,pid,43,para);
#endif

    return returnVal;
}

int   MPI_Reduce_scatter( sendbuf, recvbuf, recvcnts, datatype, op, comm ) //44
    void * sendbuf;
    void * recvbuf;
    int * recvcnts;
    MPI_Datatype datatype;
    MPI_Op op;
    MPI_Comm comm;
{
    int   returnVal;

#ifdef PERF_ASSERT
    _timer_start(PATN);
#endif 

    returnVal = PMPI_Reduce_scatter( sendbuf, recvbuf, recvcnts, datatype, op, comm );

#ifdef PERF_ASSERT
    _timer_stop(PATN);
#endif 

#ifdef PERF_ASSERT
    double r = _timer_read(PATN);
    _timer_clear(PATN);
    double e = E_MPI_Reduce_scatter( sendbuf, recvbuf, recvcnts, datatype, op, comm );
    int pid;
    PMPI_Comm_rank(MPI_COMM_WORLD, &pid);
    char para[100];
    sprintf(para, "recvcnts:%d",*recvcnts);
    R_log(R_Level,0,r,e,pid,44,para);
#endif

    return returnVal;
}

int   MPI_Reduce( sendbuf, recvbuf, count, datatype, op, root, comm ) //45
    void * sendbuf;
    void * recvbuf;
    int count;
    MPI_Datatype datatype;
    MPI_Op op;
    int root;
    MPI_Comm comm;
{
    int   returnVal;

#ifdef PERF_ASSERT
    _timer_start(PATN);
#endif 

    returnVal = PMPI_Reduce( sendbuf, recvbuf, count, datatype, op, root, comm );

#ifdef PERF_ASSERT
    _timer_stop(PATN);
#endif 

#ifdef PERF_ASSERT
    double r = _timer_read(PATN);
    _timer_clear(PATN);
    double e = E_MPI_Reduce( sendbuf, recvbuf, count, datatype, op, root, comm );
    int pid;
    PMPI_Comm_rank(MPI_COMM_WORLD, &pid);
    char para[100];
    sprintf(para, "count:%d, root:%d",count,root);
    R_log(R_Level,0,r,e,pid,45,para);
#endif

    return returnVal;
}

int   MPI_Scan( sendbuf, recvbuf, count, datatype, op, comm ) //46
    void * sendbuf;
    void * recvbuf;
    int count;
    MPI_Datatype datatype;
    MPI_Op op;
    MPI_Comm comm;
{
    int   returnVal;
    returnVal = PMPI_Scan( sendbuf, recvbuf, count, datatype, op, comm );
    return returnVal;
}

int   MPI_Scatter( sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm ) //47
    void * sendbuf;
    int sendcnt;
    MPI_Datatype sendtype;
    void * recvbuf;
    int recvcnt;
    MPI_Datatype recvtype;
    int root;
    MPI_Comm comm;
{
    int   returnVal;

#ifdef PERF_ASSERT
    _timer_start(PATN);
#endif 

    returnVal = PMPI_Scatter( sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm );

#ifdef PERF_ASSERT
    _timer_stop(PATN);
#endif 

#ifdef PERF_ASSERT
    double r = _timer_read(PATN);
    _timer_clear(PATN);
    double e = E_MPI_Scatter( sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm );
    int pid;
    PMPI_Comm_rank(MPI_COMM_WORLD, &pid);
    char para[100];
    sprintf(para, "sendcnt:%d, recvcnt:%d, root:%d",sendcnt,recvcnt,root);
    R_log(R_Level,0,r,e,pid,47,para);
#endif

    return returnVal;
}

int   MPI_Scatterv( sendbuf, sendcnts, displs, sendtype, recvbuf, recvcnt, recvtype, root, comm ) //48
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
    int returnVal;

#ifdef PERF_ASSERT
    _timer_start(PATN);
#endif 

    returnVal = PMPI_Scatterv( sendbuf, sendcnts, displs, sendtype, recvbuf, recvcnt, recvtype, root, comm );

#ifdef PERF_ASSERT
    _timer_stop(PATN);
#endif 

#ifdef PERF_ASSERT
    double r = _timer_read(PATN);
    _timer_clear(PATN);
    double e = E_MPI_Scatterv( sendbuf, sendcnts, displs, sendtype, recvbuf, recvcnt, recvtype, root, comm );
    int pid;
    PMPI_Comm_rank(MPI_COMM_WORLD, &pid);
    char para[100];
    sprintf(para, "sendcnts:%d, recvcnt:%d, root:%d",*sendcnts,recvcnt,root);
    R_log(R_Level,0,r,e,pid,48,para);
#endif

    return returnVal;
}

/******************************************************************
 *                                                                 *
 *               MPI Communicator Related Functions                *
 *                                                                 *
 ******************************************************************/

int   MPI_Comm_compare( comm1, comm2, result ) //49
    MPI_Comm comm1;
    MPI_Comm comm2;
    int * result;
{
    int   returnVal;
    returnVal = PMPI_Comm_compare( comm1, comm2, result );
    return returnVal;
}

int   MPI_Comm_create( comm, group, comm_out ) //50
    MPI_Comm comm;
    MPI_Group group;
    MPI_Comm * comm_out;
{
    int   returnVal;
    returnVal = PMPI_Comm_create( comm, group, comm_out );
    return returnVal;
}

int   MPI_Comm_dup( comm, comm_out ) //51
    MPI_Comm comm;
    MPI_Comm * comm_out;
{
    int   returnVal;
    returnVal = PMPI_Comm_dup( comm, comm_out );
    return returnVal;
}

int   MPI_Comm_free( comm ) //52
    MPI_Comm * comm;
{
    int   returnVal;
    returnVal = PMPI_Comm_free( comm );
    return returnVal;
}

int   MPI_Comm_group( comm, group ) //53
    MPI_Comm comm;
    MPI_Group * group;
{
    int   returnVal;
    returnVal = PMPI_Comm_group( comm, group );
    return returnVal;
}

int   MPI_Comm_rank( comm, rank ) //54
    MPI_Comm comm;
    int * rank;
{
    int   returnVal;
    returnVal = PMPI_Comm_rank( comm, rank );
    return returnVal;
}

int   MPI_Comm_remote_group( comm, group ) //55
    MPI_Comm comm;
    MPI_Group * group;
{
    int   returnVal;
    returnVal = PMPI_Comm_remote_group( comm, group );
    return returnVal;
}

int   MPI_Comm_remote_size( comm, size ) //56
    MPI_Comm comm;
    int * size;
{
    int   returnVal;
    returnVal = PMPI_Comm_remote_size( comm, size );
    return returnVal;
}

int   MPI_Comm_size( comm, size ) //57
    MPI_Comm comm;
    int * size;
{
    int   returnVal;
    returnVal = PMPI_Comm_size( comm, size );
    return returnVal;
}

int   MPI_Comm_split( comm, color, key, comm_out ) //58
    MPI_Comm comm;
    int color;
    int key;
    MPI_Comm * comm_out;
{
    int   returnVal;
    returnVal = PMPI_Comm_split( comm, color, key, comm_out );
    return returnVal;
}

int   MPI_Comm_test_inter( comm, flag ) //59
    MPI_Comm comm;
    int * flag;
{
    int   returnVal;
    returnVal = PMPI_Comm_test_inter( comm, flag );
    return returnVal;
}

int   MPI_Group_compare( group1, group2, result ) //60
    MPI_Group group1;
    MPI_Group group2;
    int * result;
{
    int   returnVal;
    returnVal = PMPI_Group_compare( group1, group2, result );
    return returnVal;
}

int   MPI_Group_difference( group1, group2, group_out ) //61
    MPI_Group group1;
    MPI_Group group2;
    MPI_Group * group_out;
{
    int   returnVal;
    returnVal = PMPI_Group_difference( group1, group2, group_out );
    return returnVal;
}

int   MPI_Group_excl( group, n, ranks, newgroup ) //62
    MPI_Group group;
    int n;
    int * ranks;
    MPI_Group * newgroup;
{
    int   returnVal;
    returnVal = PMPI_Group_excl( group, n, ranks, newgroup );
    return returnVal;
}

int   MPI_Group_free( group ) //63
    MPI_Group * group;
{
    int   returnVal;
    returnVal = PMPI_Group_free( group );
    return returnVal;
}

int   MPI_Group_incl( group, n, ranks, group_out ) //64
    MPI_Group group;
    int n;
    int * ranks;
    MPI_Group * group_out;
{
    int   returnVal;
    returnVal = PMPI_Group_incl( group, n, ranks, group_out );
    return returnVal;
}

int   MPI_Group_intersection( group1, group2, group_out ) //65
    MPI_Group group1;
    MPI_Group group2;
    MPI_Group * group_out;
{
    int   returnVal;
    returnVal = PMPI_Group_intersection( group1, group2, group_out );
    return returnVal;
}

int   MPI_Group_rank( group, rank ) //66
    MPI_Group group;
    int * rank;
{
    int   returnVal;
    returnVal = PMPI_Group_rank( group, rank );
    return returnVal;
}

int   MPI_Group_range_excl( group, n, ranges, newgroup ) //67
    MPI_Group group;
    int n;
    int ranges[][3];
    MPI_Group * newgroup;
{
    int   returnVal;
    returnVal = PMPI_Group_range_excl( group, n, ranges, newgroup );
    return returnVal;
}

int   MPI_Group_range_incl( group, n, ranges, newgroup ) //68
    MPI_Group group;
    int n;
    int ranges[][3];
    MPI_Group * newgroup;
{
    int   returnVal;
    returnVal = PMPI_Group_range_incl( group, n, ranges, newgroup );
    return returnVal;
}

int   MPI_Group_size( group, size ) //69
    MPI_Group group;
    int * size;
{
    int   returnVal;
    returnVal = PMPI_Group_size( group, size );
    return returnVal;
}

int   MPI_Group_translate_ranks( group_a, n, ranks_a, group_b, ranks_b ) //70
    MPI_Group group_a;
    int n;
    int * ranks_a;
    MPI_Group group_b;
    int * ranks_b;
{
    int   returnVal;
    returnVal = PMPI_Group_translate_ranks( group_a, n, ranks_a, group_b, ranks_b );
    return returnVal;
}

int   MPI_Group_union( group1, group2, group_out ) //71
    MPI_Group group1;
    MPI_Group group2;
    MPI_Group * group_out;
{
    int   returnVal;
    returnVal = PMPI_Group_union( group1, group2, group_out );
    return returnVal;
}

int   MPI_Intercomm_create( local_comm, local_leader, peer_comm, remote_leader, tag, comm_out ) //72
    MPI_Comm local_comm;
    int local_leader;
    MPI_Comm peer_comm;
    int remote_leader;
    int tag;
    MPI_Comm * comm_out;
{
    int   returnVal;
    returnVal = PMPI_Intercomm_create( local_comm, local_leader, peer_comm, remote_leader, tag, comm_out );
    return returnVal;
}

int   MPI_Intercomm_merge( comm, high, comm_out ) //73
    MPI_Comm comm;
    int high;
    MPI_Comm * comm_out;
{
    int   returnVal;
    returnVal = PMPI_Intercomm_merge( comm, high, comm_out );
    return returnVal;
}

int   MPI_Keyval_create( copy_fn, delete_fn, keyval, extra_state ) //74
    MPI_Copy_function * copy_fn;
    MPI_Delete_function * delete_fn;
    int * keyval;
    void * extra_state;
{
    int   returnVal;
    returnVal = PMPI_Keyval_create( copy_fn, delete_fn, keyval, extra_state );
    return returnVal;
}

int   MPI_Keyval_free( keyval ) //75
    int * keyval;
{
    int   returnVal;
    returnVal = PMPI_Keyval_free( keyval );
    return returnVal;
}

int   MPI_Cart_coords( comm, rank, maxdims, coords ) //76
    MPI_Comm comm;
    int rank;
    int maxdims;
    int * coords;
{
    int   returnVal;
    returnVal = PMPI_Cart_coords( comm, rank, maxdims, coords );
    return returnVal;
}

int   MPI_Cart_create( comm_old, ndims, dims, periods, reorder, comm_cart ) //77
    MPI_Comm comm_old;
    int ndims;
    int * dims;
    int * periods;
    int reorder;
    MPI_Comm * comm_cart;
{
    int   returnVal;
    returnVal = PMPI_Cart_create( comm_old, ndims, dims, periods, reorder, comm_cart );
    return returnVal;
}

int   MPI_Cart_get( comm, maxdims, dims, periods, coords ) //78
    MPI_Comm comm;
    int maxdims;
    int * dims;
    int * periods;
    int * coords;
{
    int   returnVal;
    returnVal = PMPI_Cart_get( comm, maxdims, dims, periods, coords );
    return returnVal;
}

int   MPI_Cart_map( comm_old, ndims, dims, periods, newrank ) //79
    MPI_Comm comm_old;
    int ndims;
    int * dims;
    int * periods;
    int * newrank;
{
    int   returnVal;
    returnVal = PMPI_Cart_map( comm_old, ndims, dims, periods, newrank );
    return returnVal;
}

int   MPI_Cart_rank( comm, coords, rank ) //80
    MPI_Comm comm;
    int * coords;
    int * rank;
{
    int   returnVal;
    returnVal = PMPI_Cart_rank( comm, coords, rank );
    return returnVal;
}

int   MPI_Cart_shift( comm, direction, displ, source, dest ) //81
    MPI_Comm comm;
    int direction;
    int displ;
    int * source;
    int * dest;
{
    int   returnVal;
    returnVal = PMPI_Cart_shift( comm, direction, displ, source, dest );
    return returnVal;
}

int   MPI_Cart_sub( comm, remain_dims, comm_new ) //82
    MPI_Comm comm;
    int * remain_dims;
    MPI_Comm * comm_new;
{
    int   returnVal;
    returnVal = PMPI_Cart_sub( comm, remain_dims, comm_new );
    return returnVal;
}

int   MPI_Cartdim_get( comm, ndims ) //83
    MPI_Comm comm;
    int * ndims;
{
    int   returnVal;
    returnVal = PMPI_Cartdim_get( comm, ndims );
    return returnVal;
}

int  MPI_Dims_create( nnodes, ndims, dims ) //84
    int nnodes;
    int ndims;
    int * dims;
{
    int  returnVal;
    returnVal = PMPI_Dims_create( nnodes, ndims, dims );
    return returnVal;
}

int   MPI_Graph_create( comm_old, nnodes, index, edges, reorder, comm_graph ) //85
    MPI_Comm comm_old;
    int nnodes;
    int * index;
    int * edges;
    int reorder;
    MPI_Comm * comm_graph;
{
    int   returnVal;
    returnVal = PMPI_Graph_create( comm_old, nnodes, index, edges, reorder, comm_graph );
    return returnVal;
}

int   MPI_Graph_get( comm, maxindex, maxedges, index, edges ) //86
    MPI_Comm comm;
    int maxindex;
    int maxedges;
    int * index;
    int * edges;
{
    int   returnVal;
    returnVal = PMPI_Graph_get( comm, maxindex, maxedges, index, edges );
    return returnVal;
}

int   MPI_Graph_map( comm_old, nnodes, index, edges, newrank ) //87
    MPI_Comm comm_old;
    int nnodes;
    int * index;
    int * edges;
    int * newrank;
{
    int   returnVal;
    returnVal = PMPI_Graph_map( comm_old, nnodes, index, edges, newrank );
    return returnVal;
}

int   MPI_Graph_neighbors( comm, rank, maxneighbors, neighbors ) //88
    MPI_Comm comm;
    int rank;
    int  maxneighbors;
    int * neighbors;
{
    int   returnVal;
    returnVal = PMPI_Graph_neighbors( comm, rank, maxneighbors, neighbors );
    return returnVal;
}

int   MPI_Graph_neighbors_count( comm, rank, nneighbors ) //89
    MPI_Comm comm;
    int rank;
    int * nneighbors;
{
    int   returnVal;
    returnVal = PMPI_Graph_neighbors_count( comm, rank, nneighbors );
    return returnVal;
}

int   MPI_Graphdims_get( comm, nnodes, nedges ) //90
    MPI_Comm comm;
    int * nnodes;
    int * nedges;
{
    int   returnVal;
    returnVal = PMPI_Graphdims_get( comm, nnodes, nedges );
    return returnVal;
}

int   MPI_Topo_test( comm, top_type ) //91
    MPI_Comm comm;
    int * top_type;
{
    int   returnVal;
    returnVal = PMPI_Topo_test( comm, top_type );
    return returnVal;
}

/******************************************************************
 *                                                                 *
 *                      MPI Other Functions                        *
 *                                                                 *
 ******************************************************************/

/* LAM MPI defines MPI_Abort as a macro! We check for this and if 
   it is defined that way, we change the MPI_Abort wrapper */

#if (defined(MPI_Abort) && defined(_ULM_MPI_H_))
int _MPI_Abort( MPI_Comm comm, int errorcode, char * file, int line) //92
#else
int  MPI_Abort( comm, errorcode )
    MPI_Comm comm;
    int errorcode;
#endif /* MPI_Abort & LAM MPI [LAM MPI] */
{
    int  returnVal;
    returnVal = PMPI_Abort( comm, errorcode );
    return returnVal;
}

int  MPI_Error_class( errorcode, errorclass ) //93
    int errorcode;
    int * errorclass;
{
    int  returnVal;
    returnVal = PMPI_Error_class( errorcode, errorclass );
    return returnVal;
}

int  MPI_Errhandler_create( function, errhandler ) //94
    MPI_Handler_function * function;
    MPI_Errhandler * errhandler;
{
    int  returnVal;
    returnVal = PMPI_Errhandler_create( function, errhandler );
    return returnVal;
}

int  MPI_Errhandler_free( errhandler ) //95
    MPI_Errhandler * errhandler;
{
    int  returnVal;
    returnVal = PMPI_Errhandler_free( errhandler );
    return returnVal;
}

int  MPI_Errhandler_get( comm, errhandler ) //96
    MPI_Comm comm;
    MPI_Errhandler * errhandler;
{
    int  returnVal;
    returnVal = PMPI_Errhandler_get( comm, errhandler );
    return returnVal;
}

int  MPI_Error_string( errorcode, string, resultlen ) //97
    int errorcode;
    char * string;
    int * resultlen;
{
    int  returnVal;
    returnVal = PMPI_Error_string( errorcode, string, resultlen );
    return returnVal;
}

int  MPI_Errhandler_set( comm, errhandler ) //98
    MPI_Comm comm;
    MPI_Errhandler errhandler;
{
    int  returnVal;
    returnVal = PMPI_Errhandler_set( comm, errhandler );
    return returnVal;
}


int  MPI_Get_processor_name( name, resultlen ) //99
    char * name;
    int * resultlen;
{
    int  returnVal;
    returnVal = PMPI_Get_processor_name( name, resultlen );
    return returnVal;
}


double  MPI_Wtick(  ) //100
{
    double  returnVal;
    returnVal = PMPI_Wtick(  );
    return returnVal;
}

double  MPI_Wtime(  ) //101
{
    double  returnVal;
    returnVal = PMPI_Wtime(  );
    return returnVal;
}

int  MPI_Address( location, address ) //102
    void * location;
    MPI_Aint * address;
{
    int  returnVal;
    returnVal = PMPI_Address( location, address );
    return returnVal;
}

int  MPI_Op_create( function, commute, op ) //103
    MPI_User_function * function;
    int commute;
    MPI_Op * op;
{
    int  returnVal;
    returnVal = PMPI_Op_create( function, commute, op );
    return returnVal;
}

int  MPI_Op_free( op ) //104
    MPI_Op * op;
{
    int  returnVal;
    returnVal = PMPI_Op_free( op );
    return returnVal;
}

int   MPI_Attr_delete( comm, keyval ) //105
    MPI_Comm comm;
    int keyval;
{
    int   returnVal;
    returnVal = PMPI_Attr_delete( comm, keyval );
    return returnVal;
}

int   MPI_Attr_get( comm, keyval, attr_value, flag ) //106
    MPI_Comm comm;
    int keyval;
    void * attr_value;
    int * flag;
{
    int   returnVal;
    returnVal = PMPI_Attr_get( comm, keyval, attr_value, flag );
    return returnVal;
}

int   MPI_Attr_put( comm, keyval, attr_value ) //107
    MPI_Comm comm;
    int keyval;
    void * attr_value;
{
    int   returnVal;
    returnVal = PMPI_Attr_put( comm, keyval, attr_value );
    return returnVal;
}

int  MPI_Buffer_attach( buffer, size ) //108
    void * buffer;
    int size;
{
    int  returnVal;
    returnVal = PMPI_Buffer_attach( buffer, size );
    return returnVal;
}

int  MPI_Buffer_detach( buffer, size ) //109
    void * buffer;
    int * size;
{
    int  returnVal;
    returnVal = PMPI_Buffer_detach( buffer, size );
    return returnVal;
}

int   MPI_Get_elements( status, datatype, elements ) //110
    MPI_Status * status;
    MPI_Datatype datatype;
    int * elements;
{
    int   returnVal;
    returnVal = PMPI_Get_elements( status, datatype, elements );
    return returnVal;
}

int  MPI_Get_count( status, datatype, count ) //111
    MPI_Status * status;
    MPI_Datatype datatype;
    int * count;
{
    int  returnVal;
    returnVal = PMPI_Get_count( status, datatype, count );
    return returnVal;
}

int   MPI_Type_commit( datatype ) //112
    MPI_Datatype * datatype;
{
    int   returnVal;
    returnVal = PMPI_Type_commit( datatype );
    return returnVal;
}

int  MPI_Type_contiguous( count, old_type, newtype ) //113
    int count;
    MPI_Datatype old_type;
    MPI_Datatype * newtype;
{
    int  returnVal;
    returnVal = PMPI_Type_contiguous( count, old_type, newtype );
    return returnVal;
}

int  MPI_Type_extent( datatype, extent ) //114
    MPI_Datatype datatype;
    MPI_Aint * extent;
{
    int  returnVal;
    returnVal = PMPI_Type_extent( datatype, extent );
    return returnVal;
}

int   MPI_Type_free( datatype ) //115
    MPI_Datatype * datatype;
{
    int   returnVal;
    returnVal = PMPI_Type_free( datatype );
    return returnVal;
}

int  MPI_Type_hindexed( count, blocklens, indices, old_type, newtype ) //116
    int count;
    int * blocklens;
    MPI_Aint * indices;
    MPI_Datatype old_type;
    MPI_Datatype * newtype;
{
    int  returnVal;
    returnVal = PMPI_Type_hindexed( count, blocklens, indices, old_type, newtype );
    return returnVal;
}

int  MPI_Type_hvector( count, blocklen, stride, old_type, newtype ) //117
    int count;
    int blocklen;
    MPI_Aint stride;
    MPI_Datatype old_type;
    MPI_Datatype * newtype;
{
    int  returnVal;
    returnVal = PMPI_Type_hvector( count, blocklen, stride, old_type, newtype );
    return returnVal;
}

int  MPI_Type_indexed( count, blocklens, indices, old_type, newtype ) //118
    int count;
    int * blocklens;
    int * indices;
    MPI_Datatype old_type;
    MPI_Datatype * newtype;
{
    int  returnVal;
    returnVal = PMPI_Type_indexed( count, blocklens, indices, old_type, newtype );
    return returnVal;
}

int   MPI_Type_lb( datatype, displacement ) //119
    MPI_Datatype datatype;
    MPI_Aint * displacement;
{
    int   returnVal;
    returnVal = PMPI_Type_lb( datatype, displacement );
    return returnVal;
}

int   MPI_Type_size( datatype, size ) //120
    MPI_Datatype datatype;
    int * size;
{
    int   returnVal;
    returnVal = PMPI_Type_size( datatype, size );
    return returnVal;
}

int  MPI_Type_struct( count, blocklens, indices, old_types, newtype ) //121
    int count;
    int * blocklens;
    MPI_Aint * indices;
    MPI_Datatype * old_types;
    MPI_Datatype * newtype;
{
    int  returnVal;
    returnVal = PMPI_Type_struct( count, blocklens, indices, old_types, newtype );
    return returnVal;
}

int   MPI_Type_ub( datatype, displacement ) //122
    MPI_Datatype datatype;
    MPI_Aint * displacement;
{
    int   returnVal;
    returnVal = PMPI_Type_ub( datatype, displacement );
    return returnVal;
}

int  MPI_Type_vector( count, blocklen, stride, old_type, newtype ) //123
    int count;
    int blocklen;
    int stride;
    MPI_Datatype old_type;
    MPI_Datatype * newtype;
{
    int  returnVal;
    returnVal = PMPI_Type_vector( count, blocklen, stride, old_type, newtype );
    return returnVal;
}

int   MPI_Unpack( inbuf, insize, position, outbuf, outcount, type, comm ) //124
    void * inbuf;
    int insize;
    int * position;
    void * outbuf;
    int outcount;
    MPI_Datatype type;
    MPI_Comm comm;
{
    int   returnVal;
    returnVal = PMPI_Unpack( inbuf, insize, position, outbuf, outcount, type, comm );
    return returnVal;
}

int   MPI_Pack( inbuf, incount, type, outbuf, outcount, position, comm ) //125
    void * inbuf;
    int incount;
    MPI_Datatype type;
    void * outbuf;
    int outcount;
    int * position;
    MPI_Comm comm;
{
    int   returnVal;
    returnVal = PMPI_Pack( inbuf, incount, type, outbuf, outcount, position, comm );
    return returnVal;
}

int   MPI_Pack_size( incount, datatype, comm, size ) //126
    int incount;
    MPI_Datatype datatype;
    MPI_Comm comm;
    int * size;
{
    int   returnVal;
    returnVal = PMPI_Pack_size( incount, datatype, comm, size );
    return returnVal;
}


