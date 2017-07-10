#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_ADAPTIVEGREEDY_H
#define LAWA_METHODS_ADAPTIVE_SOLVERS_ADAPTIVEGREEDY_H 1

#include <vector>

#include <engine.h>

#include <lawa/methods/adaptive/datastructures/htcoefficients.h>
#include <lawa/methods/adaptive/operators/operatorsd/sepop.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/solvers/solver_parameters.h>
#include <lawa/righthandsides/separablerhsd.h>

namespace lawa
{

template <typename Optype, typename Basis, typename Prec, typename T>
unsigned
adaptiveGreedy(      Engine*                          ep,
                     Sepop<Optype>&                   A,
                     Sepdiagscal<Basis>&              S,
                     Prec&                            P,
                     HTCoefficients<T, Basis>&        u,
                     std::vector<IndexSet<Index1D> >& Lambda,
                     SeparableRHSD<T, Basis>&         f,
               const AdaptiveGreedyParams&            params);


} // namespace lawa

#include <lawa/methods/adaptive/solvers/adaptiveGreedy.tcc>

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_ADAPTIVEGREEDY_H
