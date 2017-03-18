#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_RANK1ALS_H
#define LAWA_METHODS_ADAPTIVE_SOLVERS_RANK1ALS_H 1

#include <vector>

#include <lawa/methods/adaptive/datastructures/htcoefficients.h>
#include <lawa/methods/adaptive/operators/operatorsd/sepop.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/preconditioners/preconditionersd/sepdiagscal.h>

namespace lawa
{

/* Rank 1 one symmetric (A s.p.d.) solver over fixed wavelet index set
 * Requires valid initial x
 * Uses cg for leaf optimization
 */
template <typename Optype, typename T, typename Basis>
unsigned
rank1als_sym(       Sepop<Optype>&                      A,
                    HTCoefficients<T, Basis>&           x,
              const HTCoefficients<T, Basis>&           b,
              const std::vector<IndexSet<Index1D> >&    Lambda,
                    T&                                  residual,
              const bool                                orthog    = false,
              const T                                   tol       = 1e-08,
              const unsigned                            max_sweep = 100,
              const T                                   tol_cg    = 1e-08,
              const unsigned                            maxit_cg  = 1e+06);


/* Preconditioned version */
template <typename Optype, typename T, typename Basis>
unsigned
rank1als_sym(       Sepop<Optype>&                      A,
                    Sepdiagscal<Basis>&                 S,
                    HTCoefficients<T, Basis>&           x,
              const HTCoefficients<T, Basis>&           b,
              const std::vector<IndexSet<Index1D> >&    Lambda,
                    T&                                  residual,
              const bool                                orthog    = false,
              const T                                   tol       = 1e-08,
              const unsigned                            max_sweep = 100,
              const T                                   tol_cg    = 1e-08,
              const unsigned                            maxit_cg  = 1e+06);


/* Simple preconditioning */
template <typename Optype, typename Prec, typename T, typename Basis>
unsigned
rank1als_sym(       Sepop<Optype>&                      A,
                    Prec&                               P,
                    HTCoefficients<T, Basis>&           x,
              const HTCoefficients<T, Basis>&           b,
              const std::vector<IndexSet<Index1D> >&    Lambda,
                    T&                                  residual,
              const bool                                orthog    = false,
              const T                                   tol       = 1e-08,
              const unsigned                            max_sweep = 100,
              const T                                   tol_cg    = 1e-08,
              const unsigned                            maxit_cg  = 1e+06);


/* Rank 1 preconditioning */
template <typename Optype, typename Prec, typename T, typename Basis>
unsigned
precrank1als_sym(       Sepop<Optype>&                      A,
                        Prec&                               P,
                        HTCoefficients<T, Basis>&           x,
                  const HTCoefficients<T, Basis>&           b,
                  const std::vector<IndexSet<Index1D> >&    Lambda,
                        T&                                  residual,
                  const bool                                check_res = false,
                  const bool                                orthog    = false,
                  const bool                                sw        = true,
                  const T                                   balance   = 1.,
                  const T                                   tol       = 1e-08,
                  const unsigned                            max_sweep = 100,
                  const T                                   tol_cg    = 1e-08,
                  const unsigned                            maxit_cg  = 1e+06);

} // namespace lawa

#include <lawa/methods/adaptive/solvers/rank1als.tcc>

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_RANK1ALS_H
