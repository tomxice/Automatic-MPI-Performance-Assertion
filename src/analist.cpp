#include<iostream>
#include<map>
#include<vector>
#include<string>
#include<cstdio>
#include<cstdlib>
#include<cstring>

using namespace std;

class FuncRecord {
public:
    string name;
    string location;
    string para;
    double t_exp;
    int sc, nc, lc;
    FuncRecord() {
        sc = nc = lc = 0;
    }
    bool operator = (const FuncRecord &r) {
        name = r.name;
        location = r.location;
        para = r.para;
        t_exp = r.t_exp;
        sc = r.sc, nc = r.nc, lc = r.lc;
    }
    bool operator < (const FuncRecord &r) const {
        if (name != r.name) return name < r.name;
        if (location != r.location) return location < r.location;
        if (t_exp != r.t_exp) return t_exp < r.t_exp;
        return false;
    }
    bool operator == (const FuncRecord &r) const {
        return name == r.name && location == r.location && t_exp == r.t_exp;
    }
    bool operator != (const FuncRecord &r) const {
        return ! (*this == r);
    }
};

const char* call_begin = "Func_name";
map<FuncRecord, vector<double> > info;
FILE* pfinal;

void parse_file(int file_num, int i) {
    char fn[100];
    sprintf(fn, "mpipa.report.%d", i);
    FILE* pf = fopen(fn,"r");
    if (!pf) {
        cout << "file open failed!" << endl;
    }
    char line[1000];
    char nnn[128];
    while (fgets(line, 1000, pf)) {
        sscanf(line, "%s", nnn);
        if (string(nnn) == call_begin)
            break;
    }
    //printf("Begin !!\n");
    FuncRecord rec;
    rec.name = "null";
    FuncRecord tmp;
    while (true) {
        char* cont = fgets(line, 1000, pf);
        if (cont == NULL) {
            if (rec.name != "null") {
                if (info.find(rec) == info.end()) {
                    info[rec].resize(file_num);
                    for (int k = 0; k < file_num; ++ k)
                        info[rec][k] = -1;
                }
                //printf("lc:%d, sc:%d, nc:%d \n",rec.lc, rec.sc, rec.nc);
                info[rec][i] = ((double)rec.lc+rec.sc)/(rec.lc+rec.sc+rec.nc);
                //printf("%s:%d:%lf\n",rec.name.c_str(),i,info[rec][i]);
            }
            break;
        }
        char func[28], la[40], lb[40], lc[40], type[20],para[128];
        double t_exp, nouse;
        int count;
        sscanf(line, "%s%s%s%s%s%lf%d%lf%lf%lf%s", func, la, lb, lc, type, &t_exp, &count, &nouse, &nouse, &nouse, para);  
        tmp.name = string(func);
        tmp.location = string(la);
        tmp.location.append(" ");
        tmp.location.append(lb);
        tmp.location.append(" ");
        tmp.location.append(lc);
        tmp.t_exp = t_exp;
        tmp.para = para;
        if (tmp != rec) {
            if (rec.name != "null") {
                if (info.find(rec) == info.end()) {
                    info[rec].resize(file_num);
                    for (int k = 0; k < file_num; ++ k)
                        info[rec][k] = -1;
                }
                //printf("lc:%d, sc:%d, nc:%d \n",rec.lc, rec.sc, rec.nc);
                info[rec][i] = ((double)rec.lc+rec.sc)/(rec.lc+rec.sc+rec.nc);
                //printf("%s:%d:%lf\n",rec.name.c_str(),i,info[rec][i]);
            }
            rec = tmp; 
        }
        if (strcmp(type,"NORMAL") == 0) rec.nc = count;
        if (strcmp(type,"SHORT") == 0) rec.sc = count;
        if (strcmp(type,"LONG") == 0) rec.lc = count;
    }
    fclose(pf);
}

//decendant
void sort_p(vector<double> &a){
    int size = a.size();
    for (int i = 0; i < size; ++ i) {
        double m = i;
        for (int j = i; j < size; ++ j) {
            if (a[j] > a[m]) m = j;
        }
        double tmp = a[m];
        a[m] = a[i];
        a[i] = tmp;
    }
}

void analyse(map<FuncRecord, vector<double> >::iterator iter) {
    vector<double> &a = iter->second;
    sort_p(a);
#if 0
    for (int i = 0; i < a.size(); ++ i) {
        fprintf(pfinal, "%lf ",a[i]);
    }
    fprintf(pfinal, "\n");
#endif
    double p_avg, max_avg, min_avg;
    double th_ideal = 1e-5;
    double th_noise = 0.1;
    double th_lb = 2;
    double sum = 0;
    int size;
    int length = a.size();
    for (size = 0; size < length; ++ size) {
        if (a[size] < 0) break;  
        sum += a[size];
    }
    p_avg = sum/size;
    if (p_avg < th_ideal) {
        fprintf(pfinal,"%-28s%-40s%-20s%-s\n",iter->first.name.c_str(),iter->first.location.c_str(),"GOOD JOB",iter->first.para.c_str());
        return;
    }
    if (p_avg < th_noise) {
        fprintf(pfinal,"%-28s%-40s%-20s%-s\n",iter->first.name.c_str(),iter->first.location.c_str(),"SYSTEM NOISE",iter->first.para.c_str());
        return;
    }
    double gap=-1, tmp;
    int index = -1;
    for (int i = 0; i < size-1; ++ i) {
        tmp = a[i] - a[i+1];
        if (tmp > gap) {
            gap = tmp;
            index = i;
        }
    }
    index ++;
    max_avg = 0;
    for (int i = 0; i < index; ++ i) {
        max_avg += a[i];
    }
    max_avg /= index;
    min_avg = 0;
    for (int i = index; i < size; ++ i) {
        min_avg += a[i];
    }
    min_avg /= (size-index);
    if (max_avg/min_avg > th_lb) {
        fprintf(pfinal,"%-28s%-40s%-20s%-s\n",iter->first.name.c_str(),iter->first.location.c_str(),"LOAD IMBALANCE",iter->first.para.c_str());
    }
    else {
        fprintf(pfinal,"%-28s%-40s%-20s%-s\n",iter->first.name.c_str(),iter->first.location.c_str(),"CONTENTION",iter->first.para.c_str());
        //fprintf(pfinal,"avg:%.4lf,max_avg:%.4lf,min_avg:%.4lf\n",p_avg,max_avg,min_avg);
    }
}
    


int main(int argc, char** argv) {
    if (argc != 2) {
        cout << "Usage: analist.out procs" << endl;
        return -1;
    }
    int file_num = atoi(argv[1]);
    for (int i = 0; i < file_num; ++ i) {
        parse_file(file_num, i);
    }
    pfinal = fopen("final.txt","w+");
    map<FuncRecord, vector<double> >::iterator iter;
    for (iter = info.begin(); iter != info.end(); ++ iter) {
        analyse(iter);
    }
    fclose(pfinal);
    return 0;
}
