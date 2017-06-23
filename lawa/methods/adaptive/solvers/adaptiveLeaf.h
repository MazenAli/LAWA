#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_adaptiveLeaf_H
#define LAWA_METHODS_ADAPTIVE_SOLVERS_adaptiveLeaf_H 1

#include <flens/flens.cxx>

#include <lawa/methods/adaptive/operators/operatorsd/sepop.h>
#include <lawa/methods/adaptive/datastructures/htcoefficients.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/index.h>

namespace lawa
{

// This assumes a pre-specified structure of the equation!
// See thesis.
// Assumes rhs is given in CP form and b is already computed over
// Lambda in CP form.
template <typename Optype, typename T, typename Prec, typename Rhs>
unsigned
adaptiveLeaf_laplace(      Sepop<Optype>&                              A,
                     const flens::GeMatrix<
                           flens::FullStorage< T, flens::ColMajor > >& Pj,
                           Prec&                                       P,
                           HTCoefficients<T, Basis>&                   u,
                           IndexSet<Index1D>&                          Lambda,
                     const SeparableRHSD<T, Basis>&                    f,
                     const HTCoefficients<T, Basis>&                   b,
                     const unsigned                                    j,
                     const AdaptiveAlsParams&                          p);

} // namespace lawa

#include <lawa/methods/adaptive/solvers/adaptiveLeaf.tcc>

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_adaptiveLeaf_H
