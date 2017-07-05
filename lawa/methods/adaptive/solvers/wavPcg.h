#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_WAVPCG_H
#define LAWA_METHODS_ADAPTIVE_SOLVERS_WAVPCG_H 1

#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/datastructures/index.h>

namespace lawa
{

// Pcg solver for the ALS routines
template <typename Operator, typename Preconditioner, typename T>
unsigned
wavPcg(      Operator&                                  A,
             Preconditioner&                            P,
             Coefficients<Lexicographical, T, Index1D>& x,
       const Coefficients<Lexicographical, T, Index1D>& b,
             T&                                         res_cg,
       const T                                          tol     = 1e-08,
       const unsigned                                   maxit   = 1e+05,
       const bool                                       verbose = false);

} // namespace lawa

#include <lawa/methods/adaptive/solvers/wavPcg.tcc>
#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_WAVPCG_H
