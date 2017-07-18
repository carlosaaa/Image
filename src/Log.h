#ifndef __LOG_H__
#define __LOG_H__

#include <stdio.h>

extern bool g_log_open;
extern FILE* g_fp;

#define LOG(FMT, ...) if (g_log_open) fprintf(g_fp, FMT, __VA_ARGS__)



#endif /* __LOG_H__ */

