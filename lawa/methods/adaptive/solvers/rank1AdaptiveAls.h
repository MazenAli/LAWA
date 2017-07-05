#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_RANK1ADAPTIVEALS_H
#define LAWA_METHODS_ADAPTIVE_SOLVERS_RANK1ADAPTIVEALS_H 1

#include <vector>

#include <lawa/methods/adaptive/operators/operatorsd/sepop.h>
#include <lawa/methods/adaptive/datastructures/htcoefficients.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/solvers/solver_parameters.h>

namespace lawa
{

// Computes a rank 1 solution via an ALS procedure with adaptive discretization
template <typename Optype, typename Prec,
          typename T, typename Basis, typename Rhs>
unsigned
rank1AdaptiveAls(      Sepop<Optype>&                                A,
                       Prec&                                         P,
                       SepCoefficients<Lexicographical, T, Index1D>& v,
                       HTCoefficients<T, Basis>&                     x,
                       Rhs&                                          f,
                       std::vector<IndexSet<Index1D> >&              Lambda,
                 const Rank1AdaptiveAlsParams&                       params);


} // namespace lawa

#include <lawa/methods/adaptive/solvers/rank1AdaptiveAls.tcc>

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_RANK1ADAPTIVEALS_H
