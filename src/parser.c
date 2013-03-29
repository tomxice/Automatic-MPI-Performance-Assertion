#include <stdlib.h>
#include <stdio.h>

#define VERIFY 0

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
typedef struct IMBPara {
    int n_barrier;
    int n_bcast;
    int n_reduce;
    int n_gather;
    int n_allreduce;
    int n_allgather;
    struct Para {
        int bytes;
        double t_min, t_max, t_avg;
    };
    struct BarrierPara {
        int proc; 
        int n_byte;
        double t_min, t_max, t_avg;
    } barrier[PROC_N];
    struct BcastPara {
        int proc;
        int n_byte;
        struct Para para[BYTE_N];
    } bcast[PROC_N];
    struct ReducePara {
        int proc;
        int n_byte;
        struct Para para[BYTE_N];
    } reduce[PROC_N];
    struct GatherPara {
        int proc;
        int n_byte;
        struct Para para[BYTE_N];
    } gather[PROC_N];
    struct AllreducePara {
        int proc;
        int n_byte;
        struct Para para[BYTE_N];
    } allreduce[PROC_N];
    struct AllgatherPara {
        int proc;
        int n_byte;
        struct Para para[BYTE_N];
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
    data,
    unknown
} IMBState;


// loggpo measure points --> MPS
// maybe 1000+, which size varys
// from 0 to 1,000,000+, size may
// not use as index.
#define MPS 1500
LoggpoPara logps[MPS]; 

IMBPara imb;

// parameters
// f_para: the data file
// loggps: an array of struct where to store data,
//         pass the head pointer here
void parse_loggpo(const char* f_para, pLoggpoPara logps) {
    FILE* pf_para = NULL;
    pf_para = fopen(f_para,"r");
    if (! pf_para) {
        printf("file open failed!\n");
    }

    // read the file
    int index = 0;
    int LINES = 200;
    char line[LINES];
    while (line == fgets(line, LINES, pf_para)) {
        pLoggpoPara p = &(logps[index]); 
        int ret = sscanf(line, "%d%lf%lf%lf%lf%lf%lf%lf", 
          &p->size, &p->os, &p->or, &p->ov, &p->sr, &p->gap, &p->rtt, &p->rtt100);
        if (ret == 8) ++index;
    }
    // verify
    #if VERIFY
    for (int i = 0; i < index; ++ i) {
        pLoggpoPara p = &(logps[i]); 
        printf("\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
          p->size, p->os, p->or, p->ov, p->sr, p->gap, p->rtt, p->rtt100);
    }
    #endif
}

IMBState test_benchmark(const char* line) {
    int true = 1, false = 0;
    int comment = false;
    for (int i = 0; i < strlen(line); ++ i) {
        if (line[i] == '#') {
            comment = true;
            //printf("this line is commented\n");
            break;
        }
    }
    if (comment == false) return data;

    char str0[100], str1[100], str2[100];
    int ret = sscanf(line, "%s%s%s", str0, str1, str2);
    //printf("ret=%d, str0:%s\n",ret,str0);
    if (ret != 3) return init;
    //printf ("sscanf ret == 3\n");
    //printf ("|%s| |%s| |%s| \n", str0,str1,str2);
    if (strcmp(str1,"Benchmarking") != 0) return init;

    IMBState s = unknown;
    if (strcmp(str2,"Barrier") == 0) s = barrier;
    else if (strcmp(str2,"Bcast") == 0) s = bcast;
    else if (strcmp(str2,"Reduce") == 0) s = reduce;
    else if (strcmp(str2,"Gather") == 0) s = gather;
    else if (strcmp(str2,"Allreduce") == 0) s = allreduce;
    else if (strcmp(str2,"Allgather") == 0) s = allgather;
    else s = unknown;
    return s;
}
    

void parse_imb(const char* f_para, pIMBPara pimb) {
    memset(pimb, 0, sizeof(IMBPara));
    FILE* pf_para = NULL;
    pf_para = fopen(f_para,"r");
    if (! pf_para) {
        printf("file open failed!\n");
    }

    // read the file
    int proc=-1;
    int LINES = 200;
    char line[LINES];
    char str0[100], str1[100], str2[100];
    IMBState s = init, ns;
    // read a line
    while (line == fgets(line, LINES, pf_para)) {
        printf("%s,%d\n",line,s);
        // empty line
        if (strlen(line) == 0) {
            ns = data;
            continue;
        }
            
        ns = test_benchmark(line);
        switch (s) {
        case init:
            if (ns >= init && ns <= allgather)
                s = ns;
            break;
        case barrier:
            if (ns == init) {
                int ret = sscanf(line, "%s%s%s%d", str0, str1, str2, &proc);
                if (ret != 4) break;
                if (strcmp(str1,"#processes") == 0)
                    pimb->barrier[pimb->n_barrier].proc = proc;
            }
            else if (ns == data) {
                int rep;
                double t0, t1, t2;
                int ret = sscanf(line,"%d%lf%lf%lf", &rep, &t0, &t1, &t2);
                if (ret != 4) {
                    //printf("line:%s,ret=%d, %d %f %f %f\n", line,ret, rep, t0, t1, t2);
                    //printf("Read parameters error\n");
                    break;
                }
                pimb->barrier[pimb->n_barrier].t_min = t0;
                pimb->barrier[pimb->n_barrier].t_max = t1;
                pimb->barrier[pimb->n_barrier].t_avg = t2;
            }
            else if (ns == barrier) {
                ++ pimb->n_barrier;
            }
            else {
                s = ns;
            }
            break;
        case bcast:
            if (ns == init) {
                int ret = sscanf(line, "%s%s%s%d", str0, str1, str2, &proc);
                if (ret != 4) break;
                if (strcmp(str1,"#processes") == 0)
                    pimb->bcast[pimb->n_bcast].proc = proc;
            }
            else if (ns == data) {
                int bytes, rep;
                double t0, t1, t2;
                int ret = sscanf(line,"%d%d%lf%lf%lf", &bytes, &rep, &t0, &t1, &t2);
                if (ret != 5) {
                    printf("Read parameters error\n");
                    break;
                }
                int n_byte = pimb->bcast[pimb->n_bcast].n_byte;
                pimb->bcast[pimb->n_bcast].para[n_byte].bytes = bytes;
                pimb->bcast[pimb->n_bcast].para[n_byte].t_min = t0;
                pimb->bcast[pimb->n_bcast].para[n_byte].t_max = t1;
                pimb->bcast[pimb->n_bcast].para[n_byte].t_avg = t2;
            }
            else if (ns == bcast ) {
                ++ pimb->n_bcast;
            }
            else {
                s = ns;
            }
            break;
        case reduce:
            break;
        case gather:
            break;
        case allreduce:
            break;
        case allgather:
            break;
        case data:
        case unknown:
        default:
            printf("Entering impossible State!!\n");
            break;
        }
    }
}

int main() {
    //parse_loggpo("net_para", logps);
    parse_imb("coll_para", &imb);
    return 0;
}

          
        
