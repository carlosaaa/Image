#include "Log.h"


FILE* __init()
{
	return fopen("log.txt", "w+");
}


FILE* g_fp = __init();
bool g_log_open = true;
