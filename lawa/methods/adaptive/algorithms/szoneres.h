#ifndef LAWA_METHODS_ADAPTIVE_ALGORITHMS_SZONERES_H
#define LAWA_METHODS_ADAPTIVE_ALGORITHMS_SZONERES_H 1

#include <vector>

#include <lawa/methods/adaptive/datastructures/htcoefficients.h>
#include <lawa/methods/adaptive/datastructures/sepcoefficients.h>
#include <lawa/methods/adaptive/operators/operatorsd/sepop.h>
#include <lawa/righthandsides/separablerhsd.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>

namespace lawa
{

/* Approximate residual evaluation via security zone in each dimension */
template <typename Optype, typename Basis, typename T>
HTCoefficients<T, Basis>
szoneres(      Sepop<Optype>& A,
               HTCoefficients<T, Basis>& u,
               HTCoefficients<T, Basis>& f,
               SepCoefficients<Lexicographical, T, Index1D>& fcp,
         const SeparableRHSD<T, Basis>& fint,
         const std::vector<IndexSet<Index1D> >& current,
               std::vector<IndexSet<Index1D> >& sweep,
               std::vector<IndexSet<Index1D> >& total);

} // namespace lawa

#include <lawa/methods/adaptive/algorithms/szoneres.tcc>

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_SZONERES_H
