#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_HTAWGM_H
#define LAWA_METHODS_ADAPTIVE_SOLVERS_HTAWGM_H 1

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

/* Adaptive accuracy scheme implemented as in Tobler algorithm 9 */
template <typename Optype, typename Basis, typename T>
unsigned
galerkin_pcg(Sepop<Optype>& A,
             Sepdiagscal<Basis>& S,
             HTCoefficients<T, Basis>& x,
             const HTCoefficients<T, Basis>& b,
             const std::vector<IndexSet<Index1D> >& Lambda,
             T& residual,
             const bool uzero     = true,
             const T tol          = 1e-08,
             const unsigned maxit = 1e+02,
             const T delta1       = 1e-01,
             const T delta2       = 1e-01,
             const T delta3       = 1e-01,
             const T trunc        = 1e-10);


/* The version with S^{-1}AS^{-1} */
template <typename Optype, typename Basis, typename T>
unsigned
galerkin_pcg2(Sepop<Optype>& A,
              Sepdiagscal<Basis>& S,
              HTCoefficients<T, Basis>& x,
              const HTCoefficients<T, Basis>& b,
              const std::vector<IndexSet<Index1D> >& Lambda,
              T& residual,
              const bool uzero     = true,
              const T tol          = 1e-08,
              const unsigned maxit = 1e+02,
              const T delta        = 1e-01,
              const T dres         = 2.,
              const T trunc        = 1e-10);


template <typename Optype, typename Basis, typename T>
std::vector<IndexSet<Index1D> >
presidual(Sepop<Optype>& A,
          Sepdiagscal<Basis>& S,
                HTCoefficients<T, Basis>& u,
                HTCoefficients<T, Basis>& f,
                SepCoefficients<Lexicographical, T, Index1D>& fcp,
                HTCoefficients<T, Basis>& r,
          const SeparableRHSD<T, Basis>& fint,
          const std::vector<IndexSet<Index1D> >& current,
          const std::vector<IndexSet<Index1D> >& sweep,
                std::vector<IndexSet<Index1D> >& total,
          const T trunc      = 1e-08);


/* The version with S^{-1}AS^{-1} */
template <typename Optype, typename Basis, typename T>
std::vector<IndexSet<Index1D> >
presidual2(Sepop<Optype>& A,
           Sepdiagscal<Basis>& S,
                 HTCoefficients<T, Basis>& u,
                 HTCoefficients<T, Basis>& f,
                 SepCoefficients<Lexicographical, T, Index1D>& fcp,
                 HTCoefficients<T, Basis>& r,
           const SeparableRHSD<T, Basis>& fint,
           const std::vector<IndexSet<Index1D> >& current,
           const std::vector<IndexSet<Index1D> >& sweep,
                 std::vector<IndexSet<Index1D> >& total,
           const T trunc      = 1e-08);


template <typename T, typename Basis>
std::vector<IndexSet<Index1D> >
bulk(const T alpha, const T resex,
           HTCoefficients<T, Basis>& res,
           std::vector<IndexSet<Index1D> >& Lambda,
     const std::vector<IndexSet<Index1D> >& diff);


template <typename T, typename Basis>
std::vector<IndexSet<Index1D> >
bulkBestN(const T alpha, const T resex,
           HTCoefficients<T, Basis>& res,
           std::vector<IndexSet<Index1D> >& Lambda,
     const std::vector<IndexSet<Index1D> >& diff);


template <typename Optype, typename Basis, typename T>
unsigned
htawgm(Sepop<Optype>&                   A,
       Sepdiagscal<Basis>&              S,
       HTCoefficients<T, Basis>&        u,
       const SeparableRHSD<T, Basis>&   f,
       std::vector<IndexSet<Index1D> >& Lambda,
       double&                          residual,
       HTAWGM_Params&                   params);


/* The version with S^{-1}AS^{-1} */
template <typename Optype, typename Basis, typename T>
unsigned
htawgm2(      Sepop<Optype>&                   A,
              Sepdiagscal<Basis>&              S,
              HTCoefficients<T, Basis>&        u,
        const SeparableRHSD<T, Basis>&   f,
              std::vector<IndexSet<Index1D> >& Lambda,
              double&                          residual,
              HTAWGM_Params&                   params);

} // namespace lawa

#include <lawa/methods/adaptive/solvers/htawgm.tcc>

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_HTAWGM_H
