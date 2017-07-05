#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_adaptiveLeaf_H
#define LAWA_METHODS_ADAPTIVE_SOLVERS_adaptiveLeaf_H 1

#include <vector>

#include <flens/flens.cxx>

#include <lawa/methods/adaptive/operators/operatorsd/sepop.h>
#include <lawa/methods/adaptive/datastructures/htcoefficients.h>
#include <lawa/methods/adaptive/datastructures/sepcoefficients.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/solvers/solver_parameters.h>
#include <lawa/righthandsides/separablerhsd.h>

namespace lawa
{

// This assumes a pre-specified structure of the equation!
// See thesis.
// Assumes rhs is given in CP form and b is already computed over
// Lambda in CP form.
template <typename Optype, typename Prec, typename T, typename Basis,
          typename Rhs>
unsigned
adaptiveLeaf(      Sepop<Optype>&                                A,
                   Prec&                                         P,
                   HTCoefficients<T, Basis>&                     u,
                   Coefficients<Lexicographical, T, Index1D>&    v,
                   std::vector<IndexSet<Index1D>>&               active,
                   Rhs&                                          f,
             const unsigned                                      j,
             const AdaptiveLeafParams&                           params);

} // namespace lawa

#include <lawa/methods/adaptive/solvers/adaptiveLeaf.tcc>

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_adaptiveLeaf_H
