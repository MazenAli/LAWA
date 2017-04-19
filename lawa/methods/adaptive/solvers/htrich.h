#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_HTRICH_H
#define LAWA_METHODS_ADAPTIVE_SOLVERS_HTRICH_H 1

#include <vector>

#include <lawa/methods/adaptive/datastructures/htcoefficients.h>
#include <lawa/methods/adaptive/operators/operatorsd/sepop.h>
#include <lawa/methods/adaptive/preconditioners/preconditionersd/sepdiagscal.h>
#include <lawa/righthandsides/separablerhsd.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/solvers/solver_parameters.h>

namespace lawa
{

/* Implementation according to Dahmen/Bachmayr */
template <typename Optype, typename Basis, typename T>
unsigned
htrich(      Sepop<Optype>&                   A,
             Sepdiagscal<Basis>&              S,
             HTCoefficients<T, Basis>&        u,
       const SeparableRHSD<T, Basis>&         f,
             std::vector<IndexSet<Index1D> >& Lambda,
             T&                               residual,
       const HTRICH_Params&                   params);

template <typename T, typename Basis>
std::vector<IndexSet<Index1D> >
coarsen(      HTCoefficients<T, Basis>&         u,
        const std::vector<IndexSet<Index1D> >&  Lambda,
        const T                                 eps);

template <typename T>
unsigned
findminj(const T rho, const T omega, const T beta1,
         const T beta2, const T kappa1,
         unsigned max = 100);

} // namespace lawa

#include <lawa/methods/adaptive/solvers/htrich.tcc>

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_HTRICH_H
