#ifndef LAWA_EXTENSIONS_MATLAB_INTERFACE_CPP
#define LAWA_EXTENSIONS_MATLAB_INTERFACE_CPP 1

#include <cassert>
#include <cstring>

#include "interface.h"

namespace flens
{
namespace matlab
{

mxArray*
create2DCell(const std::vector<std::vector<RealColMat>>& B)
{
    assert(B.size()>0);
    assert(B[0].size()>0);

    typedef RealColMat::ElementType ElementType;

    auto k = B.size();
    auto d = B[0].size();
    mxArray *icell;
    icell = mxCreateCellMatrix(d, k);
    mwIndex ij[2];

    for (mwIndex i=0; i<d; ++i) {
        for (mwIndex j=0; j<k; ++j) {
            assert(B[j].size()==d);

            auto m     = B[j][i].numRows();
            auto n     = B[j][i].numCols();
            mxArray *A = mxCreateDoubleMatrix(m, n, mxREAL);
            std::memcpy((void*) mxGetPr(A), (void*) B[j][i].data(),
                        m*n*sizeof(ElementType));

            ij[0] = i;
            ij[1] = j;
            mwIndex index = mxCalcSingleSubscript(icell, k, ij);
            mxSetCell(icell, index, A);
        }
    }

    return icell;
}

std::vector<std::vector<RealColMat>>
read2DCell(const mxArray* cell)
{
    typedef RealColMat::ElementType ElementType;

    const mwSize *dimCell = mxGetDimensions(cell);
    const mwSize numdim   = mxGetNumberOfDimensions(cell);
    assert(numdim==2);
    (void) numdim;

    auto d = dimCell[0];
    auto k = dimCell[1];
    std::vector<std::vector<RealColMat>>  ret(k);
    mwIndex ij[2];
    for (mwIndex j=0; j<k; ++j) {
        ret[j].resize(d);
        for (mwIndex i=0; i<d; ++i) {
            ij[0]         = i;
            ij[1]         = j;
            mwIndex index = mxCalcSingleSubscript(cell, k, ij);
            mxArray *A    = mxGetCell(cell, index);

            const mwSize *dimA = mxGetDimensions(A);
            ret[j][i].resize(dimA[0], dimA[1]);
            std::memcpy((void*) ret[j][i].data(), (void*) mxGetPr(A),
                        dimA[0]*dimA[1]*sizeof(ElementType));
        }
    }

    return ret;
}

} // namespace matlab
} // namespace flens

#endif // LAWA_EXTENSIONS_MATLAB_INTERFACE_CPP
