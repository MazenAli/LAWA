#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_SPLITTING_H
#define LAWA_METHODS_ADAPTIVE_SOLVERS_SPLITTING_H 1

#include <vector>
#include <flens/flens.cxx>

#include <lawa/methods/adaptive/operators/operatorsd/sepop.h>
#include <lawa/methods/adaptive/datastructures/htcoefficients.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>

namespace lawa
{

template <typename Optype, typename Prec, typename T, typename Basis>
unsigned
symm_splitting(  Sepop<Optype>&               A,
                 Prec&                        P,
                 flens::GeMatrix
                 <flens::FullStorage
                 <T, cxxblas::ColMajor> >&    V,
           const flens::GeMatrix
                   <flens::FullStorage
                 <T, cxxblas::ColMajor> >&    B,
           const flens::GeMatrix
                 <flens::FullStorage
                 <T, cxxblas::ColMajor> >&    Pj,
                 HTCoefficients<T, Basis>&    u,
           const unsigned                     j,
           const IndexSet<Index1D>&           Lambda,
                 T&                           res,
           const T                            tol    = 1e-08,
           const unsigned                     maxit  = 1e+06);

} // namespace lawa

#include <lawa/methods/adaptive/solvers/splitting.tcc>

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_SPLITTING_H
