#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "mpiname.h"

#define PATN 1

FILE* f_detail;

int E_result[sizeof(MPI_Functions)][3];
//TODO variable array in C ?

void E_report(double r, double e, int id);

/******************************************************************
*                                                                 *
*                   MPI Functions for Management                  *
*                                                                 *
******************************************************************/
double E_MPI_Init(
int * argc,
char *** argv );

#ifdef PERF_MPI_THREADED
double E_MPI_Init_thread (
int * argc,
char *** argv,
int required,
int *provided );

#endif /* PERF_MPI_THREADED */
double E_MPI_Finalize();


/******************************************************************
*                                                                 *
*          MPI Point-to-Point(
*                                                                 *
******************************************************************/
double E_MPI_Bsend(
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int tag,
MPI_Comm comm );

double E_MPI_Bsend_init(
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int tag,
MPI_Comm comm,
MPI_Request * request );

double E_MPI_Recv_init(
void * buf,
int count,
MPI_Datatype datatype,
int source,
int tag,
MPI_Comm comm,
MPI_Request * request );

double E_MPI_Send_init(
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int tag,
MPI_Comm comm,
MPI_Request * request );

double E_MPI_Ibsend(
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int tag,
MPI_Comm comm,
MPI_Request * request );
double E_MPI_Irecv(
void * buf,
int count,
MPI_Datatype datatype,
int source,
int tag,
MPI_Comm comm,
MPI_Request * request );
double E_MPI_Irsend(
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int tag,
MPI_Comm comm,
MPI_Request * request );
double E_MPI_Isend(
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int tag,
MPI_Comm comm,
MPI_Request * request );
double E_MPI_Issend(
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int tag,
MPI_Comm comm,
MPI_Request * request );
double E_MPI_Recv(
void * buf,
int count,
MPI_Datatype datatype,
int source,
int tag,
MPI_Comm comm,
MPI_Status * status );
double E_MPI_Rsend(
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int tag,
MPI_Comm comm );
double E_MPI_Rsend_init(
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int tag,
MPI_Comm comm,
MPI_Request * request );
double E_MPI_Send(
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int tag,
MPI_Comm comm );
double E_MPI_Sendrecv(
void * sendbuf,
int sendcount,
MPI_Datatype sendtype,
int dest,
int sendtag,
void * recvbuf,
int recvcount,
MPI_Datatype recvtype,
int source,
int recvtag,
MPI_Comm comm,
MPI_Status * status );
double E_MPI_Sendrecv_replace(
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int sendtag,
int source,
int recvtag,
MPI_Comm comm,
MPI_Status * status );
double E_MPI_Ssend(
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int tag,
MPI_Comm comm );
double E_MPI_Ssend_init(
void * buf,
int count,
MPI_Datatype datatype,
int dest,
int tag,
MPI_Comm comm,
MPI_Request * request );
double E_MPI_Test(
MPI_Request * request,
int * flag,
MPI_Status * status );
double E_MPI_Testall(
int count,
MPI_Request * array_of_requests,
int * flag,
MPI_Status * array_of_statuses );
double E_MPI_Testany(
int count,
MPI_Request * array_of_requests,
int * index,
int * flag,
MPI_Status * status );
double E_MPI_Test_cancelled(
MPI_Status * status,
int * flag );
double E_MPI_Testsome(
int incount,
MPI_Request * array_of_requests,
int * outcount,
int * array_of_indices,
MPI_Status * array_of_statuses );
double E_MPI_Wait(
MPI_Request * request,
MPI_Status * status );
double E_MPI_Waitall(
int count,
MPI_Request * array_of_requests,
MPI_Status * array_of_statuses );
double E_MPI_Waitany(
int count,
MPI_Request * array_of_requests,
int * index,
MPI_Status * status );
double E_MPI_Waitsome(
int incount,
MPI_Request * array_of_requests,
int * outcount,
int * array_of_indices,
MPI_Status * array_of_statuses );
double E_MPI_Cancel(
MPI_Request * request );
double E_MPI_Request_free(
MPI_Request * request );
double E_MPI_Start(
MPI_Request * request );
double E_MPI_Startall(
int count,
MPI_Request * array_of_requests );
double E_MPI_Iprobe(
int source,
int tag,
MPI_Comm comm,
int * flag,
MPI_Status * status );
double E_MPI_Probe(
int source,
int tag,
MPI_Comm comm,
MPI_Status * status );
/******************************************************************
*                                                                 *
*               MPI Collective Communication Related              *
*                                                                 *
******************************************************************/
double E_MPI_Allgather(
void * sendbuf,
int sendcount,
MPI_Datatype sendtype,
void * recvbuf,
int recvcount,
MPI_Datatype recvtype,
MPI_Comm comm );
double E_MPI_Allgatherv(
void * sendbuf,
int sendcount,
MPI_Datatype sendtype,
void * recvbuf,
int * recvcounts,
int * displs,
MPI_Datatype recvtype,
MPI_Comm comm );
double E_MPI_Allreduce(
void * sendbuf,
void * recvbuf,
int count,
MPI_Datatype datatype,
MPI_Op op,
MPI_Comm comm );
double E_MPI_Alltoall(
void * sendbuf,
int sendcount,
MPI_Datatype sendtype,
void * recvbuf,
int recvcnt,
MPI_Datatype recvtype,
MPI_Comm comm );
double E_MPI_Alltoallv(
void * sendbuf,
int * sendcnts,
int * sdispls,
MPI_Datatype sendtype,
void * recvbuf,
int * recvcnts,
int * rdispls,
MPI_Datatype recvtype,
MPI_Comm comm );
double E_MPI_Barrier(
MPI_Comm comm );
double E_MPI_Bcast(
void * buffer,
int count,
MPI_Datatype datatype,
int root,
MPI_Comm comm );
double E_MPI_Gather(
void * sendbuf,
int sendcnt,
MPI_Datatype sendtype,
void * recvbuf,
int recvcount,
MPI_Datatype recvtype,
int root,
MPI_Comm comm );
double E_MPI_Gatherv(
void * sendbuf,
int sendcnt,
MPI_Datatype sendtype,
void * recvbuf,
int * recvcnts,
int * displs,
MPI_Datatype recvtype,
int root,
MPI_Comm comm );
double E_MPI_Reduce_scatter(
void * sendbuf,
void * recvbuf,
int * recvcnts,
MPI_Datatype datatype,
MPI_Op op,
MPI_Comm comm );
double E_MPI_Reduce(
void * sendbuf,
void * recvbuf,
int count,
MPI_Datatype datatype,
MPI_Op op,
int root,
MPI_Comm comm );
double E_MPI_Scan(
void * sendbuf,
void * recvbuf,
int count,
MPI_Datatype datatype,
MPI_Op op,
MPI_Comm comm );
double E_MPI_Scatter(
void * sendbuf,
int sendcnt,
MPI_Datatype sendtype,
void * recvbuf,
int recvcnt,
MPI_Datatype recvtype,
int root,
MPI_Comm comm );
double E_MPI_Scatterv(
void * sendbuf,
int * sendcnts,
int * displs,
MPI_Datatype sendtype,
void * recvbuf,
int recvcnt,
MPI_Datatype recvtype,
int root,
MPI_Comm comm );
/******************************************************************
*                                                                 *
*               MPI Communicator Related Functions                *
*                                                                 *
******************************************************************/
double E_MPI_Comm_compare(
MPI_Comm comm1,
MPI_Comm comm2,
int * result );
double E_MPI_Comm_create(
MPI_Comm comm,
MPI_Group group,
MPI_Comm * comm_out );
double E_MPI_Comm_dup(
MPI_Comm comm,
MPI_Comm * comm_out );
double E_MPI_Comm_free(
MPI_Comm * comm );
double E_MPI_Comm_group(
MPI_Comm comm,
MPI_Group * group );
double E_MPI_Comm_rank(
MPI_Comm comm,
int * rank );
double E_MPI_Comm_remote_group(
MPI_Comm comm,
MPI_Group * group );
double E_MPI_Comm_remote_size(
MPI_Comm comm,
int * size );
double E_MPI_Comm_size(
MPI_Comm comm,
int * size );
double E_MPI_Comm_split(
MPI_Comm comm,
int color,
int key,
MPI_Comm * comm_out );
double E_MPI_Comm_test_inter(
MPI_Comm comm,
int * flag );
double E_MPI_Group_compare(
MPI_Group group1,
MPI_Group group2,
int * result );
double E_MPI_Group_difference(
MPI_Group group1,
MPI_Group group2,
MPI_Group * group_out );
double E_MPI_Group_excl(
MPI_Group group,
int n,
int * ranks,
MPI_Group * newgroup );
double E_MPI_Group_free(
MPI_Group * group );
double E_MPI_Group_incl(
MPI_Group group,
int n,
int * ranks,
MPI_Group * group_out );
double E_MPI_Group_intersection(
MPI_Group group1,
MPI_Group group2,
MPI_Group * group_out );
double E_MPI_Group_rank(
MPI_Group group,
int * rank );
double E_MPI_Group_range_excl(
MPI_Group group,
int n,
int ranges[][3],
MPI_Group * newgroup );
double E_MPI_Group_range_incl(
MPI_Group group,
int n,
int ranges[][3],
MPI_Group * newgroup );
double E_MPI_Group_size(
MPI_Group group,
int * size );
double E_MPI_Group_translate_ranks(
MPI_Group group_a,
int n,
int * ranks_a,
MPI_Group group_b,
int * ranks_b );
double E_MPI_Group_union(
MPI_Group group1,
MPI_Group group2,
MPI_Group * group_out );
double E_MPI_Intercomm_create(
MPI_Comm local_comm,
int local_leader,
MPI_Comm peer_comm,
int remote_leader,
int tag,
MPI_Comm * comm_out );
double E_MPI_Intercomm_merge(
MPI_Comm comm,
int high,
MPI_Comm * comm_out );
double E_MPI_Keyval_create(
MPI_Copy_function * copy_fn,
MPI_Delete_function * delete_fn,
int * keyval,
void * extra_state );
double E_MPI_Keyval_free(
int * keyval );
double E_MPI_Cart_coords(
MPI_Comm comm,
int rank,
int maxdims,
int * coords );
double E_MPI_Cart_create(
MPI_Comm comm_old,
int ndims,
int * dims,
int * periods,
int reorder,
MPI_Comm * comm_cart );
double E_MPI_Cart_get(
MPI_Comm comm,
int maxdims,
int * dims,
int * periods,
int * coords );
double E_MPI_Cart_map(
MPI_Comm comm_old,
int ndims,
int * dims,
int * periods,
int * newrank );
double E_MPI_Cart_rank(
MPI_Comm comm,
int * coords,
int * rank );
double E_MPI_Cart_shift(
MPI_Comm comm,
int direction,
int displ,
int * source,
int * dest );
double E_MPI_Cart_sub(
MPI_Comm comm,
int * remain_dims,
MPI_Comm * comm_new );
double E_MPI_Cartdim_get(
MPI_Comm comm,
int * ndims );
double E_MPI_Dims_create(
int nnodes,
int ndims,
int * dims );
double E_MPI_Graph_create(
MPI_Comm comm_old,
int nnodes,
int * index,
int * edges,
int reorder,
MPI_Comm * comm_graph );
double E_MPI_Graph_get(
MPI_Comm comm,
int maxindex,
int maxedges,
int * index,
int * edges );
double E_MPI_Graph_map(
MPI_Comm comm_old,
int nnodes,
int * index,
int * edges,
int * newrank );
double E_MPI_Graph_neighbors(
MPI_Comm comm,
int rank,
int  maxneighbors,
int * neighbors );
double E_MPI_Graph_neighbors_count(
MPI_Comm comm,
int rank,
int * nneighbors );
double E_MPI_Graphdims_get(
MPI_Comm comm,
int * nnodes,
int * nedges );
double E_MPI_Topo_test(
MPI_Comm comm,
int * top_type );
/******************************************************************
*                                                                 *
*                      MPI Other Functions                        *
*                                                                 *
******************************************************************/
#if (defined(MPI_Abort) && defined(_ULM_MPI_H_))
int E_MPI_Abort( MPI_Comm comm, int errorcode, char * file, int line);
#else
int E_MPI_Abort( MPI_Comm comm, int errorcode );
#endif /* MPI_Abort & LAM MPI [LAM MPI] */


double E_MPI_Error_class(
int errorcode,
int * errorclass );
double E_MPI_Errhandler_create(
MPI_Handler_function * function,
MPI_Errhandler * errhandler );
double E_MPI_Errhandler_free(
MPI_Errhandler * errhandler );
double E_MPI_Errhandler_get(
MPI_Comm comm,
MPI_Errhandler * errhandler );
double E_MPI_Error_string(
int errorcode,
char * string,
int * resultlen );
double E_MPI_Errhandler_set(
MPI_Comm comm,
MPI_Errhandler errhandler );
double E_MPI_Get_processor_name(
char * name,
int * resultlen );
double  E_MPI_Wtick();
double  E_MPI_Wtime();
double E_MPI_Address(
void * location,
MPI_Aint * address );
double E_MPI_Op_create(
MPI_User_function * function,
int commute,
MPI_Op * op );
double E_MPI_Op_free(
MPI_Op * op );
double E_MPI_Attr_delete(
MPI_Comm comm,
int keyval );
double E_MPI_Attr_get(
MPI_Comm comm,
int keyval,
void * attr_value,
int * flag );
double E_MPI_Attr_put(
MPI_Comm comm,
int keyval,
void * attr_value );
double E_MPI_Buffer_attach(
void * buffer,
int size );
double E_MPI_Buffer_detach(
void * buffer,
int * size );
double E_MPI_Get_elements(
MPI_Status * status,
MPI_Datatype datatype,
int * elements );
double E_MPI_Get_count(
MPI_Status * status,
MPI_Datatype datatype,
int * count );
double E_MPI_Type_commit(
MPI_Datatype * datatype );
double E_MPI_Type_contiguous(
int count,
MPI_Datatype old_type,
MPI_Datatype * newtype );
double E_MPI_Type_extent(
MPI_Datatype datatype,
MPI_Aint * extent );
double E_MPI_Type_free(
MPI_Datatype * datatype );
double E_MPI_Type_hindexed(
int count,
int * blocklens,
MPI_Aint * indices,
MPI_Datatype old_type,
MPI_Datatype * newtype );
double E_MPI_Type_hvector(
int count,
int blocklen,
MPI_Aint stride,
MPI_Datatype old_type,
MPI_Datatype * newtype );
double E_MPI_Type_indexed(
int count,
int * blocklens,
int * indices,
MPI_Datatype old_type,
MPI_Datatype * newtype );
double E_MPI_Type_lb(
MPI_Datatype datatype,
MPI_Aint * displacement );
double E_MPI_Type_size(
MPI_Datatype datatype,
int * size );
double E_MPI_Type_struct(
int count,
int * blocklens,
MPI_Aint * indices,
MPI_Datatype * old_types,
MPI_Datatype * newtype );
double E_MPI_Type_ub(
MPI_Datatype datatype,
MPI_Aint * displacement );
double E_MPI_Type_vector(
int count,
int blocklen,
int stride,
MPI_Datatype old_type,
MPI_Datatype * newtype );
double E_MPI_Unpack(
void * inbuf,
int insize,
int * position,
void * outbuf,
int outcount,
MPI_Datatype type,
MPI_Comm comm );
double E_MPI_Pack(
void * inbuf,
int incount,
MPI_Datatype type,
void * outbuf,
int outcount,
int * position,
MPI_Comm comm );
double E_MPI_Pack_size(
int incount,
MPI_Datatype datatype,
MPI_Comm comm,
int * size );
/********************************************************
*                  User defined functions               *
********************************************************/
void E_MPI_Profile_on();


void E_MPI_Profile_off();


