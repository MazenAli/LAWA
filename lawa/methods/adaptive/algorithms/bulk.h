#ifndef LAWA_METHODS_ADAPTIVE_ALGORITHMS_BULK_H
#define LAWA_METHODS_ADAPTIVE_ALGORITHMS_BULK_H 1

#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>

namespace lawa
{

template <typename T, typename Basis>
IndexSet<Index1D>
bulk(      IndexSet<Index1D>&                         active,
           Coefficients<Lexicographical, T, Index1D>& r,
     const T                                          alpha,
     const T                                          res,
     const T                                          res_base,
     const Basis&                                     basis,
     const bool                                       verbose = false);

} // namespace lawa

#include <lawa/methods/adaptive/algorithms/bulk.tcc>

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_BULK_H
