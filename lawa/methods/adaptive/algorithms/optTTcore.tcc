#ifndef LAWA_METHODS_ADAPTIVE_ALGORITHMS_OPTTTCORE_TCC
#define LAWA_METHODS_ADAPTIVE_ALGORITHMS_OPTTTCORE_TCC 1

#include <iostream>
#include <cassert>

#include <extensions/matlab/matlab.h>

namespace lawa
{

template <typename T, typename I>
std::vector<flens::GeMatrix<
flens::FullStorage<T, cxxblas::ColMajor> > >
optTTcoreLaplace(      Engine                                       *ep,
                 const std::vector<flens::GeMatrix<
                       flens::FullStorage<T, cxxblas::ColMajor> > >& A,
                 const std::vector<flens::GeMatrix<
                       flens::FullStorage<T, cxxblas::ColMajor> > >& rhs,
                 const std::vector<flens::GeMatrix<
                       flens::FullStorage<T, cxxblas::ColMajor> > >& x0,
                       flens::DenseVector<
                       flens::Array<I> >&                            ranks,
                 const OptTTCoreParams&                              params)
{
    assert(ep);
    assert(A.size()>0);
    assert(rhs.size()>0);
    assert(x0.size()>0);
    assert(A.size()==rhs.size());
    assert(x0.size()==A.size()-1);
    assert((unsigned) ranks.length()==x0.size());

    typedef typename flens::
                     GeMatrix
                     <flens::FullStorage<T, cxxblas::ColMajor> > Matrix;

    #ifdef VERBOSE
        char buffer[BUFSIZE];
    #endif

    mxArray *maxIt    = flens::matlab::createDouble((double) params.maxIt);
    mxArray *tol      = flens::matlab::createDouble((double) params.tol);
    mxArray *stag     = flens::matlab::createDouble((double) params.stag);
    mxArray *cellA    = flens::matlab::create1DCell(A);
    mxArray *cellrhs  = flens::matlab::create1DCell(rhs);
    mxArray *cellx0   = flens::matlab::create1DCell(x0);
    mxArray *vecranks = flens::matlab::createIndexVector(ranks);

    /* Pass input */
    (void) engPutVariable(ep, "Alawa", cellA);
    (void) engPutVariable(ep, "blawa", cellrhs);
    (void) engPutVariable(ep, "x0lawa", cellx0);
    (void) engPutVariable(ep, "r0lawa", vecranks);
    (void) engPutVariable(ep, "maxIt", maxIt);
    (void) engPutVariable(ep, "tol", tol);
    (void) engPutVariable(ep, "stag", stag);

    /* Call function */
    #ifdef VERBOSE
        (void) engOutputBuffer(ep, buffer, BUFSIZE);
    #endif
    (void)
    matlab::doEval
    (ep,
     "addpath(genpath(getenv('MATLAB_APPROXIMATIONTB')))");
    #ifdef VERBOSE
        (void) std::printf(buffer);
    #endif
    (void)
    matlab::doEval
    (ep,
    "[xlawa, rlawa] = LAWAoptTTcoreLaplace(Alawa, blawa, x0lawa, r0lawa, maxIt, tol, stag);");
    #ifdef VERBOSE
        (void) std::printf(buffer);
    #endif

    /* Process result */
    mxArray *xlawa = nullptr;
    mxArray *rlawa = nullptr;
    xlawa          = engGetVariable(ep, "xlawa");
    rlawa          = engGetVariable(ep, "rlawa");
    if (!xlawa || !rlawa) {
        std::cerr <<
        "optTTcoreLaplace: Reading variables from workspace failed\n";
    }

    std::vector<Matrix> ret    = flens::matlab::read1DCell(xlawa);
                        ranks  = flens::matlab::readIndexVector(rlawa);

    /* Clean up */
    mxDestroyArray(cellA);
    mxDestroyArray(cellrhs);
    mxDestroyArray(cellx0);
    mxDestroyArray(xlawa);
    mxDestroyArray(rlawa);

    return ret;
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_OPTTTCORE_TCC
