//
// Created by danielsaad on 4/11/17.
//
#include <fstream>
#include "util.hpp"


#ifdef MYDEBUG
ofstream dbg_file("dbg.log");
#endif

#ifdef REPORT
ofstream report_file("report.log");
#endif

#ifdef MEM_MONITOR
mem_monitor mm("mem-mon-out.csv");
#endif