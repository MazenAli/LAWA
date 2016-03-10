#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_INDEXOPS_H
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_INDEXOPS_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>

namespace lawa
{

template <typename Index, typename Basis>
unsigned long
maptoint(const Index&, const Basis&);


template <typename Basis>
unsigned long
maptoint(const Index1D& index, const Basis& basis);


template <typename Index, typename Basis>
unsigned long
maxintind(const IndexSet<Index>& indexset, const Basis& basis);


template <typename Index, typename Basis>
Index
maptowav(const unsigned long, const Basis&);


template <typename Basis>
Index1D
maptowav(const unsigned long i, const Basis& basis);

} // namespace lawa

#include <lawa/methods/adaptive/datastructures/indexops.tcc>

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_INDEXOPS_H
