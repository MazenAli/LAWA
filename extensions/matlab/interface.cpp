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
createDouble(const double d)
{
    mxArray *out = mxCreateDoubleMatrix(1, 1, mxREAL);
    double  *p   = mxGetPr(out);
    *p = d;

    return out;
}

mxArray*
create1DCell(const std::vector<RealColMat>& B)
{
    assert(B.size()>0);

    typedef RealColMat::ElementType ElementType;

    auto d = B.size();
    mxArray *icell;
    icell = mxCreateCellMatrix(d, 1);

    mwIndex id[2];
    for (mwIndex i=0; i<d; ++i) {
        auto m     = B[i].numRows();
        auto n     = B[i].numCols();
        mxArray *A = mxCreateDoubleMatrix(m, n, mxREAL);
        std::memcpy((void*) mxGetPr(A), (void*) B[i].data(),
                    m*n*sizeof(ElementType));

        id[0] = i;
        id[1] = 1;
        mwIndex index = mxCalcSingleSubscript(icell, 1, id);
        mxSetCell(icell, index, A);
    }

    return icell;
}


mxArray*
createIndexVector(const IndexVector& x)
{
    assert(x.length()>0);

    typedef flens::DenseVector<flens::Array<__UINT64_TYPE__>> ConversionType;

    auto m      = x.length();
    mwSize dims[2];
    dims[0] = m;
    dims[1] = 1;
    mxArray *A  = mxCreateNumericArray(2, dims, mxUINT64_CLASS, mxREAL);
    assert(A);
    ConversionType x_ = x;
    std::memcpy((void*) mxGetData(A), (void*) x_.data(),
                m*sizeof(__UINT64_TYPE__));
    return A;
}


std::vector<RealColMat>
read1DCell(const mxArray* cell)
{
    assert(cell);

    typedef RealColMat::ElementType ElementType;

    const mwSize *dimCell = mxGetDimensions(cell);
    const mwSize numdim   = mxGetNumberOfDimensions(cell);
    assert(numdim==2);
    assert(dimCell[1]==1);
    (void) numdim;

    auto d = dimCell[0];
    std::vector<RealColMat>  ret(d);
    mwIndex id[2];
    for (mwIndex i=0; i<d; ++i) {
        id[0]         = i;
        id[1]         = 1;
        mwIndex index = mxCalcSingleSubscript(cell, 1, id);
        mxArray *A    = mxGetCell(cell, index);

        assert(A);
        const mwSize *dimA = mxGetDimensions(A);
        ret[i].resize(dimA[0], dimA[1]);
        std::memcpy((void*) ret[i].data(), (void*) mxGetPr(A),
                    dimA[0]*dimA[1]*sizeof(ElementType));
    }

    return ret;
}

IndexVector
readIndexVector(const mxArray* vector)
{
    assert(vector);

    typedef flens::DenseVector<flens::Array<__UINT64_TYPE__>> ConversionType;

    const mwSize *dimCell = mxGetDimensions(vector);
    const mwSize numdim   = mxGetNumberOfDimensions(vector);
    assert(numdim==2);
    assert(dimCell[1]==1);
    (void) numdim;

    auto d = dimCell[0];
    IndexVector     ret(d);
    ConversionType  ret_(d);
    std::memcpy((void*) ret_.data(), (void*) mxGetData(vector),
                d*sizeof(__UINT64_TYPE__));
    ret = ret_;

    return ret;
}

} // namespace matlab
} // namespace flens

#endif // LAWA_EXTENSIONS_MATLAB_INTERFACE_CPP
