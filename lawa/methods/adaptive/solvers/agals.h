#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_AGALS_H
#define LAWA_METHODS_ADAPTIVE_SOLVERS_AGALS_H 1

#include <vector>

#include <lawa/methods/adaptive/datastructures/htcoefficients.h>
#include <lawa/methods/adaptive/operators/operatorsd/sepop.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/solvers/solver_parameters.h>

#include <engine.h>

namespace lawa
{

/* Adaptive Greedy ALS Galerkin Solver
 * So far works only with TT tree (Anthony's optTTcore)
 */
template <typename Optype, typename Prec, typename T, typename Basis>
unsigned
agals_sym(            Engine                             *ep,
                      Sepop<Optype>&                      A,
                      Prec&                               P,
                      HTCoefficients<T, Basis>&           x,
                const HTCoefficients<T, Basis>&           b,
                      std::vector<IndexSet<Index1D> >&    Lambda,
                      T&                                  residual,
                const AgALSParams&                        params);

} // namespace lawa

#include <lawa/methods/adaptive/solvers/agals.tcc>

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_AGALS_H
