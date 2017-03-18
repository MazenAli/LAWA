#ifndef LAWA_EXTENSIONS_MATLAB_DOEVAL_H
#define LAWA_EXTENSIONS_MATLAB_DOEVAL_H 1

#include <string>

#include <engine.h>

#ifndef ERRBUFSIZE
    #define ERRBUFSIZE 4096
#endif

namespace matlab
{

/* Wrapper to catch error messages */
int
doEval(Engine *ep, std::string com);

} // namespace matlab

#endif // LAWA_EXTENSIONS_MATLAB_DOEVAL_H
