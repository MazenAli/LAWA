#ifndef EXTENSIONS_GENERAL_VMEMUSAGE_CPP
#define EXTENSIONS_GENERAL_VMEMUSAGE_CPP 1

#include "vmemusage.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

namespace extensions
{

int
parseLine(char* line)
{
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

/* Note: this value is in KB! */
int
vmemusage()
{
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmSize:", 7) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}

double
kbtogb(int kb)
{
    double gb = (double) kb/1024/1024;
    return gb;
}

} // namespace extensions

#endif // EXTENSIONS_GENERAL_VMEMUSAGE_CPP
