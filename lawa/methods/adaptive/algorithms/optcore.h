#ifndef LAWA_METHODS_ADAPTIVE_ALGORITHMS_OPTCORE_H
#define LAWA_METHODS_ADAPTIVE_ALGORITHMS_OPTCORE_H 1

#include <vector>
#include <engine.h>
#include <flens/flens.cxx>

namespace lawa
{

template <typename T>
std::vector<std::vector<flens::GeMatrix<
flens::FullStorage<T, cxxblas::ColMajor> > > >
optcore(Engine *ep,
        const std::vector<std::vector<flens::GeMatrix<
              flens::FullStorage<T, cxxblas::ColMajor> > > >& B);

} // namespace lawa

#include <lawa/methods/adaptive/algorithms/optcore.tcc>

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_OPTCORE_H
