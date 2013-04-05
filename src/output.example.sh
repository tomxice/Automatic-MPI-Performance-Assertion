Automatic MPI Performance Assertion Tools Report
================================================
#There should be information about this project and copyrights

#Summary gives a brief report.
#How many call in total
#How many calls are incredible fast
#How many calls are as expected
#How many calls are stick...
Summary
==============

Function Summary
----------------
Function    Total   Short   Normal  Long
ALL         15      1       10      4
MPI_Send    6       0       4       2
MPI_Recv    6       0       6       0
MPI_Bcast   3       1       0       2

Process Summary
----------------
Process     Total   Short   Normal  Long
ALL         15      1       10      4
proc_0      9       1       4       4
proc_1      6       0       6       0

#Exception Section picks all abnormal events
#Which function is bad, and where is this functions
#What may be the reason(given in Detect field)
#Also log the arguments, expectation and real_time.
Exception
==============

Func_Name   Location    Detect          t_exp   t_real  parameter
MPI_Send    sender.c:35 Noise           4us     34us    from:proc_0 to:proc_1 ....
MPI_Send    sender.c:35 Noise           ...
MPI_Bcast   sender.c:91 Load_Balance    ...

#The Detail Section looks like a trace log
#It provide rich info to debug
Detail
==============

Function Detail
----------------
Func_Name   Location    Total   Exception
MPI_Send    sender.c:35 3       2 
**Function Trace**
No. PID t_real  t_exp   Detect  t_start t_end   parameter   
0   0   34us    4us     Noise   5968us  6002us  from:proc_0 to:proc1 ..all args passed
1   0   10us    4us     Noise   6782us  6792us  from:proc_0 to:proc1 .....
2   0   3us     4us     --      8000us  8003us  from:proc_0 to:proc1 .....

Func_Name   Location    Total   Exception
MPI_Send    sender.c:47 3       0
**Function Trace**

Func_Name   Location    Total   Exception
MPI_Recv    recver.c:24 6       0
**Function Trace**

END REPORT
===============

            
