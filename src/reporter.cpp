#include <string>
#include <cstring>
#include <cstdio>
#include <vector>
#include <map>
#include "mpiname.h"
#include "reporter.h"

using namespace std;

struct Record {
    int pid;
    int status;
    double real, expc;
    string parameter;
};
map< string, vector<Record> > result[sizeof(MPI_Functions)];

#define SHORT 0
#define NORMAL 1
#define LONG 2

void R_log(int level, int warn, double real, double expc, int pid, int mpiid, const char* c_location, const char* c_para) {
    string location(c_location, strlen(c_location));
    if (level >= 4) {
        string para(c_para, strlen(c_para));
    }
    Record record;
    record.pid = pid;
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
    int numfunc = sizeof(MPI_Functions);
    int all_sum[4];
    int func_sum[numfunc][4];//short,normal,long,all
    int proc_sum[numproc][4];//short,normal,long,all
    for (int i = 0; i < numfunc; ++ i) {
        map<string, vector<Record> >::iterator iter = result[i].begin();
        while (iter != result[i].end()) {
            for (vector<Record>::iterator it = iter->second.begin();
                    it != iter->second.end(); ++ it) {
                func_sum[i][it->status] ++;
                func_sum[i][3] ++;
                proc_sum[it->pid][it->status] ++;
                proc_sum[it->pid][3] ++;
                all_sum[it->status] ++;
                all_sum[3] ++;
            }
            ++ iter;
        }
    }
    fprintf(f_final, "Summary\n");
    fprintf(f_final, "===================\n");
    fprintf(f_final, "Function Summary\n");
    fprintf(f_final, "-------------------\n");
    fprintf(f_final, "Function\t\tTotal\tShort\tNormal\tLong\n");
    fprintf(f_final, "ALL\t\t%d\t%d\t%d\t%d\n",all_sum[3],all_sum[0],all_sum[1],all_sum[2]);
    for (int i = 0; i < numfunc; ++ i) {
        fprintf(f_final, "%s\t\t%d\t%d\t%d\t%d\n",MPI_Functions[i],func_sum[i][3],func_sum[i][0],func_sum[i][1],func_sum[i][2]);
    }
    fprintf(f_final, "Process Summary\n");
    fprintf(f_final, "-------------------\n");
    fprintf(f_final, "Process\t\tTotal\tShort\tNormal\tLong\n");
    fprintf(f_final, "ALL\t\t%d\t%d\t%d\t%d\n",all_sum[3],all_sum[0],all_sum[1],all_sum[2]);
    for (int i = 0; i < numproc; ++ i) {
        fprintf(f_final, "Proc_%d\t\t%d\t%d\t%d\t%d\n",i,proc_sum[i][3],proc_sum[i][0],proc_sum[i][1],proc_sum[i][2]);
    }
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
    fprintf(f_final, "Exception\n");
    fprintf(f_final, "===================\n");
    fprintf(f_final, "Func_name\t\tLocation\tReason\tt_exp\tt_real\n");
    for (int i = 0; i < numfunc; ++ i) {
        map<string, vector<Record> >::iterator iter = result[i].begin();
        while (iter != result[i].end()) {
            for (vector<Record>::iterator it = iter->second.begin();
                    it != iter->second.end(); ++ it) {
                if (it->status == NORMAL)
                    continue;
                fprintf(f_final, "%s\t%s\t",MPI_Functions[i],iter->first.c_str());
                if (it->status == Noise)
                    fprintf(f_final, "Noise\t");
                else if (it->status == Load_UB)
                    fprintf(f_final, "Load_UB\t");
                fprintf(f_final, "%8f\t%8f\t%s\n",it->expc,it->real,it->parameter.c_str());
            }
            ++ iter;
        }
    }
    fprintf(f_final, "Func_name\t\tLocation\tReason\tt_exp\tt_real\n");
}

