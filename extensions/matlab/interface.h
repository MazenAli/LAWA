#ifndef LAWA_EXTENSIONS_MATLAB_INTERFACE_H
#define LAWA_EXTENSIONS_MATLAB_INTERFACE_H 1

#include <vector>

#include <flens/flens.cxx>

#include <matrix.h>

namespace flens
{
namespace matlab
{

typedef flens::GeMatrix<flens::FullStorage<double, cxxblas::ColMajor>>
        RealColMat;

typedef flens::DenseVector<
        flens::Array<FLENS_DEFAULT_INDEXTYPE> >
        IndexVector;

mxArray*
createDouble(const double d);

mxArray*
create1DCell(const std::vector<RealColMat>& B);

mxArray*
createIndexVector(const IndexVector& x);

std::vector<RealColMat>
read1DCell(const mxArray* cell);

IndexVector
readIndexVector(const mxArray* vector);

} // namespace matlab
} // namespace flens

#endif // LAWA_EXTENSIONS_MATLAB_INTERFACE_H
