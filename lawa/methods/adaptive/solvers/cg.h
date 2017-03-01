#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_CG_H
#define LAWA_METHODS_ADAPTIVE_SOLVERS_CG_H 1

#include <vector>
#include <flens/flens.cxx>

#include <lawa/methods/adaptive/operators/operatorsd/sepop.h>
#include <lawa/methods/adaptive/datastructures/htcoefficients.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/preconditioners/preconditionersd/sepdiagscal.h>

namespace lawa
{

/* Computes solution of reduced equation on leaf j */
template <typename Optype, typename T, typename Basis>
unsigned
cg(      Sepop<Optype>&               A,
         flens::GeMatrix
         <flens::FullStorage
         <T, cxxblas::ColMajor> >&    U,
   const flens::GeMatrix
         <flens::FullStorage
         <T, cxxblas::ColMajor> >&    B,
   const flens::GeMatrix
         <flens::FullStorage
         <T, cxxblas::ColMajor> >&    Pj,
         HTCoefficients<T, Basis>&    u,
   const unsigned                     j,
   const IndexSet<Index1D>&           Lambda,
         T&                           res_cg,
   const T                            tol_cg    = 1e-08,
   const unsigned                     maxit_cg  = 1e+06);


/* Preconditioned version */
template <typename Optype, typename T, typename Basis>
unsigned
cg(      Sepop<Optype>&               A,
         Sepdiagscal<Basis>&          S,
         flens::GeMatrix
         <flens::FullStorage
         <T, cxxblas::ColMajor> >&    U,
   const flens::GeMatrix
         <flens::FullStorage
         <T, cxxblas::ColMajor> >&    B,
   const flens::GeMatrix
         <flens::FullStorage
         <T, cxxblas::ColMajor> >&    Pj,
         HTCoefficients<T, Basis>&    u,
   const unsigned                     j,
   const IndexSet<Index1D>&           Lambda,
         T&                           res_cg,
   const T                            tol_cg    = 1e-08,
   const unsigned                     maxit_cg  = 1e+06);


/* Simple preconditioning */
template <typename Optype, typename Prec, typename T, typename Basis>
unsigned
cg(      Sepop<Optype>&               A,
         Prec&                        P,
         flens::GeMatrix
         <flens::FullStorage
         <T, cxxblas::ColMajor> >&    U,
   const flens::GeMatrix
         <flens::FullStorage
         <T, cxxblas::ColMajor> >&    B,
   const flens::GeMatrix
         <flens::FullStorage
         <T, cxxblas::ColMajor> >&    Pj,
         HTCoefficients<T, Basis>&    u,
   const unsigned                     j,
   const IndexSet<Index1D>&           Lambda,
         T&                           res_cg,
   const T                            tol_cg    = 1e-08,
   const unsigned                     maxit_cg  = 1e+06);


/* Rank 1 preconditioning */
template <typename Optype, typename Prec, typename T, typename Basis>
unsigned
cg_rank1prec(    Sepop<Optype>&               A,
                 Prec&                        P,
                 flens::GeMatrix
                 <flens::FullStorage
                 <T, cxxblas::ColMajor> >&    U,
           const flens::GeMatrix
                 <flens::FullStorage
                 <T, cxxblas::ColMajor> >&    B,
           const flens::GeMatrix
                 <flens::FullStorage
                 <T, cxxblas::ColMajor> >&    Pj,
                 HTCoefficients<T, Basis>&    u,
           const unsigned                     j,
           const IndexSet<Index1D>&           Lambda,
                 T&                           res_cg,
           const T                            tol_cg    = 1e-08,
           const unsigned                     maxit_cg  = 1e+06);

} // namespace lawa

#include <lawa/methods/adaptive/solvers/cg.tcc>

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_CG_H
