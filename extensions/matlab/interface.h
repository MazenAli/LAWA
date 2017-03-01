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

mxArray*
create2DCell(const std::vector<std::vector<RealColMat>>& B);

std::vector<std::vector<RealColMat>>
read2DCell(const mxArray* cell);

} // namespace matlab
} // namespace flens

#endif // LAWA_EXTENSIONS_MATLAB_INTERFACE_H
