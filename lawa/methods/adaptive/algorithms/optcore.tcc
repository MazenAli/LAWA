#ifndef LAWA_METHODS_ADAPTIVE_ALGORITHMS_OPTCORE_TCC
#define LAWA_METHODS_ADAPTIVE_ALGORITHMS_OPTCORE_TCC 1

#include <iostream>
#include <cassert>

#include <extensions/matlab/matlab.h>

namespace lawa
{

template <typename T>
std::vector<std::vector<flens::GeMatrix<
flens::FullStorage<T, cxxblas::ColMajor> > > >
optcore(Engine *ep,
        const std::vector<std::vector<flens::GeMatrix<
              flens::FullStorage<T, cxxblas::ColMajor> > > >& B)
{
    assert(ep);
    assert(B.size()>0);

    typedef typename flens::
                     GeMatrix
                     <flens::FullStorage<T, cxxblas::ColMajor> > GeMat;

    mxArray *icell = flens::matlab::create2DCell(B);

    /* Pass cell */
    engPutVariable(ep, "RedSys", icell);

    /* Call function (replace) */
    std::cout << "Displaying contents from MATLAB\n";
    engEvalString(ep, "disp(RedSys)");
    engEvalString(ep, "X = RedSys");

    /* Process result (replace) */
    mxArray *ocell = nullptr;
    ocell          = engGetVariable(ep, "X");
    if (!ocell) {
        std::cerr << "optcore: Reading variable from workspace failed\n";
    }

    std::vector<std::vector<GeMat>> ret = flens::matlab::read2DCell(ocell);

    /* Clean up */
    mxDestroyArray(icell);
    mxDestroyArray(ocell);

    return ret;
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_OPTCORE_TCC
