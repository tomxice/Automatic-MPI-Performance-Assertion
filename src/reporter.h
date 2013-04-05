#define R_Level 1

FILE* f_final;
FILE* f_runtime;

void R_log(int level, int warn, double real, double expc, int pid, int mpiid, const char* c_location, const char* c_para);
void R_report(int level, int numproc);
