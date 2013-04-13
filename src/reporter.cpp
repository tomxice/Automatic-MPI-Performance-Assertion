#include <execinfo.h>
#include <string>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <map>
#include "mpiname.h"
extern "C" {
#include "reporter.h"
}

using namespace std;

struct Record {
    int pid;
    int status;
    double real, expc;
    string parameter;
    Record(): pid(-1),status(-1),real(-1),expc(-1),parameter("uninit") {}
};
const int NUMFUNC=127;
map< string, vector<Record> > result[NUMFUNC];

void R_init(const char* final_file) {
    f_final = fopen(final_file,"w");
    if (f_final == NULL) {
        printf("cannot open report file\n");
    }
}

#define SCNT 5
#define ALL 0
#define SHORT 1
#define NORMAL 2
#define LONG 3
#define NODATA 4
const char* status[SCNT] = {"ALL","SHORT","NORMAL","LONG","NODATA"};

void R_log(int level, int warn, double real, double expc, int pid, int mpiid, const char* c_para) {
    // get symbol addr
    void *buffer[SCNT];
    char **addrs;
    int nptrs = backtrace(buffer,SCNT);
    addrs = backtrace_symbols(buffer, nptrs);
    string location;
    if (addrs == NULL) {
        location = "noknown";
    }
    else {
        location.assign(addrs[2]);
        free(addrs);
    }
    // end symbol addr

    if (level >= 4) {
        string para(c_para, strlen(c_para));
    }
    Record record;
    record.pid = pid;
    if (expc < 0) record.status = NODATA;
    if (real < expc*0.25) record.status = SHORT;
    else if (real > expc*4) record.status = LONG;
    else record.status = NORMAL;

    //if (record.status == NORMAL && level <= 1) 
    //    return;

    record.real = real;
    record.expc = expc;
    if (level >= 4 || record.status != NORMAL) {
        record.parameter.assign(c_para);
    }
    result[mpiid][location].push_back(record);
    
    // TODO real time warnning
}
    
void R_report(int level, int numproc) {

    // Summary Level 0
    int numfunc = NUMFUNC;
    int all_sum[SCNT];
    int func_sum[numfunc][SCNT];//short,normal,long,all
    int proc_sum[numproc][SCNT];//short,normal,long,all
    memset(all_sum,0,SCNT*sizeof(int));
    memset(func_sum,0,numfunc*SCNT*sizeof(int));
    memset(proc_sum,0,numproc*SCNT*sizeof(int));

    for (int i = 0; i < numfunc; ++ i) {
        map<string, vector<Record> >::iterator iter = result[i].begin();
        while (iter != result[i].end()) {
            for (vector<Record>::iterator it = iter->second.begin();
                    it != iter->second.end(); ++ it) {
                func_sum[i][it->status] ++;
                func_sum[i][ALL] ++;
                proc_sum[it->pid][it->status] ++;
//                printf("pid:%d status:%d\n",it->pid,it->status);
                proc_sum[it->pid][ALL] ++;
                all_sum[it->status] ++;
                all_sum[ALL] ++;
            }
            ++ iter;
        }
    }
    fprintf(f_final, "Summary\n");
    fprintf(f_final, "===================\n\n");
    fprintf(f_final, "%-27s%-12s%-12s%-12s%-12s%-12s\n", "Function","Total","Short","Normal","Long","NoData");
    fprintf(f_final, "%-27s%-12d%-12d%-12d%-12d%-12d\n", "ALL",all_sum[ALL],all_sum[SHORT],all_sum[NORMAL],all_sum[LONG],all_sum[NODATA]);
    for (int i = 0; i < numfunc; ++ i) {
        fprintf(f_final, "%-27s%-12d%-12d%-12d%-12d%-12d\n",MPI_Functions[i],func_sum[i][ALL],func_sum[i][SHORT],func_sum[i][NORMAL],func_sum[i][LONG],func_sum[i][NODATA]);
    }
//Process info is hard to collect at runtime.. 
//because procs run in diff machines..
#if 0
    fprintf(f_final, "Process Summary\n");
    fprintf(f_final, "-------------------\n");
    fprintf(f_final, "%-15s%-12s%-12s%-12s%-12s\n", "Process","Total","Short","Normal","Long");
    fprintf(f_final, "%-15s%-12d%-12d%-12d%-12d\n", "ALL",all_sum[3],all_sum[0],all_sum[1],all_sum[2]);
    for (int i = 0; i < numproc; ++ i) {
        fprintf(f_final, "Proc_%-10d%-12d%-12d%-12d%-12d\n",i,proc_sum[i][3],proc_sum[i][0],proc_sum[i][1],proc_sum[i][2]);
    }
#endif
    // End Summary
    
    // Exceptions Level 1
    // status used to present REASON here.
    // TODO no analyse yet.. all will be Noise
    const int Noise = -1, Load_UB = -2;
    for (int i = 0; i < numfunc; ++ i) {
        map<string, vector<Record> >::iterator iter = result[i].begin();
        while (iter != result[i].end()) {
            for (vector<Record>::iterator it = iter->second.begin();
                    it != iter->second.end(); ++ it) {
                it->status = Noise;//TODO analyse here
            }
            ++ iter;
        }
    }
    fprintf(f_final, "\n\nExceptions\n");
    fprintf(f_final, "===================\n\n");
    fprintf(f_final, "%-27s%-50s%-10s%-10s%-10s%-s\n","Func_name","Location","Reason","t_exp","t_real","description");
    for (int i = 0; i < numfunc; ++ i) {
        map<string, vector<Record> >::iterator iter = result[i].begin();
        while (iter != result[i].end()) {
            for (vector<Record>::iterator it = iter->second.begin();
                    it != iter->second.end(); ++ it) {
                if (it->status == NORMAL)
                    continue;
                fprintf(f_final, "%-27s%-50s",MPI_Functions[i],iter->first.c_str());
                if (it->status == Noise)
                    fprintf(f_final, "%-10s","Noise");
                else if (it->status == Load_UB)
                    fprintf(f_final, "%-10s","Unbalance");
                fprintf(f_final, "%-10f%-10f%-s\n",it->expc,it->real,it->parameter.c_str());
            }
            ++ iter;
        }
    }
}

