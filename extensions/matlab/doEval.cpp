#ifndef LAWA_EXTENSIONS_MATLAB_DOEVAL_CPP
#define LAWA_EXTENSIONS_MATLAB_DOEVAL_CPP 1

#include <cassert>
#include <cstdio>
#include <cstring>

#include "doEval.h"

namespace matlab
{

int
doEval(Engine *ep, std::string com)
{
    assert(ep);

    char buffer[ERRBUFSIZE];
    char *_com = new char[com.length()+1];
    std::strcpy(_com, com.c_str());
    std::sprintf(buffer,
    "try,%s,catch ME,ERROR_REPORT=getReport(ME),end", _com);
    delete[] _com;
    return engEvalString(ep, buffer);
}

} // namespace matlab

#endif // LAWA_EXTENSIONS_MATLAB_DOEVAL_CPP
