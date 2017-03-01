#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_RANK1UPDATE_H
#define LAWA_METHODS_ADAPTIVE_SOLVERS_RANK1UPDATE_H 1

#include <vector>

#include <lawa/methods/adaptive/datastructures/htcoefficients.h>
#include <lawa/methods/adaptive/operators/operatorsd/sepop.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/solvers/solver_parameters.h>

namespace lawa
{

/* Rank 1 one update (A s.p.d.) over fixed wavelet index set
 */
template <typename Optype, typename Prec, typename T, typename Basis>
unsigned
rank1update_sym(      Sepop<Optype>&                      A,
                      Prec&                               P,
                      HTCoefficients<T, Basis>&           x,
                const HTCoefficients<T, Basis>&           b,
                const std::vector<IndexSet<Index1D> >&    Lambda,
                const Rank1UP_Params&                     params);

} // namespace lawa

#include <lawa/methods/adaptive/solvers/rank1update.tcc>

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_RANK1UPDATE_H
