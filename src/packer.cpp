#include <iostream>
#include <sstream>
#include <cstring>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <map>
#include <vector>

using namespace std;

#define PACK_DONE 0
#define PACK_ERROR -1

struct MPI_Call {
    char func[28];
    char addr1[20], addr2[20];
    char type[10];
    char description[50];
    double t_exp;
    string location;
    bool operator < (const MPI_Call &c) const {
        if (t_exp != c.t_exp) return t_exp < c.t_exp;
        int r = strcmp(addr2, c.addr2);
        if (r != 0) return r < 0;
        r = strcmp(type, c.type);
        if (r != 0) return r < 0;
        if (location < c.location) return true;
        r = strcmp(description, c.description);
        if (r != 0) return r < 0;
        return false;
    }
};

map<string, string> addr2line;
const char* rprefix = "mpipa.log.";
const char* wprefix = "mpipa.report.";
const char* call_begin = "Func_name";
const string wrapper_f = "mpiwrapper_f.c";
char rfname[30], wfname[30];
string exefile;

string get_line(char* addr1, char* addr2) {
    bool hit1 = false;
    map<string, string>::iterator iter;
    iter = addr2line.find(string(addr1));
    if (iter != addr2line.end()) {
        hit1 = true;
        if (iter->second != wrapper_f) 
            return iter->second;
    }
    iter = addr2line.find(string(addr2));
    if (iter != addr2line.end()) {
        return iter->second;
    }

    string func, loc1, loc2;
    FILE* pp;
    char buf[128];

    sprintf(buf, "addr2line -fse %s %s", exefile.c_str(), addr1);
    //printf("Cmd: %s\n", buf);
    pp = popen(buf,"r");
    if (pp == NULL) 
        return "Error";
    if (fgets(buf, 100, pp)) func = buf;
    if (fgets(buf, 100, pp)) loc1 = buf;
    pclose(pp);
    if (loc1.substr(0,wrapper_f.size()) != wrapper_f) {
        string ret = func.substr(0,func.size()-1)+loc1.substr(0,loc1.size()-1);
        if (! hit1)
            addr2line[string(addr1)] = ret;
        return ret;
    }

    addr2line[string(addr1)] = wrapper_f;
    sprintf(buf, "addr2line -fse %s %s", exefile.c_str(), addr2);
    //printf("Cmd: %s\n", buf);
    pp = popen(buf,"r");
    if (pp == NULL) 
        return "Error";
    if (fgets(buf, 100, pp)) func = buf;
    if (fgets(buf, 100, pp)) loc2 = buf;
    pclose(pp);
    string ret = func.substr(0,func.size()-1)+loc2.substr(0,loc2.size()-1);
    addr2line[string(addr2)] = ret;
    return ret;
}

int pack_file(int i) {
    map<MPI_Call, vector<double> > trace;
    sprintf(rfname, "%s%d", rprefix, i);
    sprintf(wfname, "%s%d", wprefix, i);
    FILE* pf = fopen(rfname, "r");
    if (! pf) 
        return PACK_ERROR;
    char line[1000];
    char tmp[128];
    while (fgets(line, 1000, pf)) {
        sscanf(line, "%s", tmp);
        if (string(tmp) == call_begin) 
            break;
    }
    //printf("Begin !!\n");
    while (fgets(line, 1000, pf)) {
        MPI_Call call;
        double t_real;
        sscanf(line, "%s%s%s%s%lf%lf%s", call.func, call.addr1, call.addr2, call.type, &call.t_exp, &t_real, call.description);
        //printf("line(parsered): %s\n", line);
        call.location = get_line(call.addr1, call.addr2);
        //printf("location: %s\n", call.location.c_str());
        trace[call].push_back(t_real); 
    }
    fclose(pf);

    //printf("Analysis\n"); 
    pf = fopen(wfname, "w");
    if (! pf) 
        return PACK_ERROR;
    fprintf(pf,"%-27s%-30s%-8s%-10s%6s  %-10s%-10s%-10s%-s\n", "Func_name", "Location", "Type", "t_exp", "count", "real_avg", "real_max", "real_min", "description");
    map<MPI_Call, vector<double> >::iterator iter;
    for (iter = trace.begin(); iter != trace.end(); ++ iter) {
        double sum = 0, max = -1, min = 1.0/0.0, avg = 0;
        for (vector<double>::iterator it = iter->second.begin();
                it != iter->second.end(); ++ it) {
            sum += *it;
            if (*it > max) max = *it;
            if (*it < min) min = *it;
        }
        int count = iter->second.size();
        avg = sum/count;
        fprintf(pf,"%-27s%-30s%-8s%-10.2e%6d  %-10.2e%-10.2e%-10.2e%-s\n", iter->first.func, iter->first.location.c_str(), iter->first.type, iter->first.t_exp, count, avg, max, min, iter->first.description);
    }
    fclose(pf);
    return PACK_DONE;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        printf("Usage: ./packer <exefile> <procs>\n");
        printf("Example: ./packer cg.D.128 128\n");
        return -1;
    }
    exefile = argv[1];
    int nprocs = atoi(argv[2]);
    for (int i = 0; i < nprocs; ++ i) {
        if (pack_file(i) != PACK_DONE) {
            printf("Error while packing No.%i file\n", i);
        }
    }
    return 0;
}

