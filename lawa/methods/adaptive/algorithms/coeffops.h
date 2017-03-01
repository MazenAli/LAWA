#ifndef LAWA_METHODS_ADAPTIVE_ALGORITHMS_COEFFOPS_H
#define LAWA_METHODS_ADAPTIVE_ALGORITHMS_COEFFOPS_H 1

#include <iostream>
#include <vector>
#include <flens/flens.cxx>
#include <lawa/settings/enum.h>
#include <lawa/righthandsides/separablerhsd.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/datastructures/sepcoefficients.h>
#include <lawa/methods/adaptive/datastructures/htcoefficients.h>
#include <lawa/methods/adaptive/operators/operatorsd/sepop.h>
#include <lawa/methods/adaptive/preconditioners/preconditionersd/sepdiagscal.h>

namespace lawa
{

template <SortingCriterion S, typename T, typename Index>
void
setCoefficients(SepCoefficients<S, T, Index>& coeffs,
                const typename SepCoefficients<S, T, Index>
                ::CoeffVec& _coeffs);


template <SortingCriterion S, typename T, typename Index>
void
setCoefficients(SepCoefficients<S, T, Index>& coeffs,
                const typename SepCoefficients<S, T, Index>
                ::size_type i,
                const typename SepCoefficients<S, T, Index>
                ::size_type j,
                const typename SepCoefficients<S, T, Index>
                ::Coeff& coeff);


template <SortingCriterion S, typename T, typename Index>
void
addCoefficients(SepCoefficients<S, T, Index>& coeffs,
                const typename SepCoefficients<S, T, Index>
                ::size_type i,
                const typename SepCoefficients<S, T, Index>
                ::size_type j,
                const typename SepCoefficients<S, T, Index>
                ::Coeff& coeff);


template <SortingCriterion S, typename T, typename Index, typename Basis>
void
genCoefficients(SepCoefficients<S, T, Index>& coeffs,
                const SeparableRHSD<T, Basis>& rhs,
                const IndexSet<Index>& indexset);


template <SortingCriterion S, typename T, typename Index, typename Basis>
void
genCoefficients(SepCoefficients<S, T, Index>& coeffs,
                const SeparableRHSD<T, Basis>& rhs,
                const std::vector<IndexSet<Index>>& indexset);


template <SortingCriterion S, typename T, typename Index, typename Basis>
void
genAddCoefficients(SepCoefficients<S, T, Index>& coeffs,
                const SeparableRHSD<T, Basis>& rhs,
                const std::vector<IndexSet<Index>>& indexset);


template <SortingCriterion S, typename T, typename Index>
std::ostream& operator<<(std::ostream& s,
                         const SepCoefficients<S, T, Index>& coeffs);


template <SortingCriterion S, typename T, typename Index, typename Basis>
unsigned FLENS_DEFAULT_INDEXTYPE
maxintind(const Coefficients<S, T, Index>& coeffs, const Basis& basis);


template <SortingCriterion S, typename T, typename Index, typename Basis>
unsigned FLENS_DEFAULT_INDEXTYPE
maxintindhash(const Coefficients<S, T, Index>& coeffs,
              const int dim,
              HTCoefficients<T, Basis>& u);


template <typename T, typename Index, typename Basis>
unsigned FLENS_DEFAULT_INDEXTYPE
maxintindhash(const IndexSet<Index>& active,
              HTCoefficients<T, Basis>& u,
              const int dim);


template <typename T, SortingCriterion S, typename Index, typename Basis>
void
set(HTCoefficients<T, Basis>& tree,
    const SepCoefficients<S, T, Index>& cp);


template <typename T, typename Basis>
void
init(HTCoefficients<T, Basis>&              tree,
     const std::vector<IndexSet<Index1D> >& activex,
     const unsigned                         rank  = 1,
     const T                                value = 1.);


template <typename T, typename Basis>
void
rndinit(HTCoefficients<T, Basis>&              tree,
        const std::vector<IndexSet<Index1D> >& activex,
        const unsigned                         rank  = 1,
        const long                             seed  = -1);


template <typename T, SortingCriterion S, typename Index, typename Basis>
void
set(HTCoefficients<T, Basis>& tree, const htucker::DimensionIndex& idx,
    const unsigned FLENS_DEFAULT_INDEXTYPE col, const Coefficients<S, T, Index>& coeff);


template <typename T, SortingCriterion S, typename Index, typename Basis>
void
set(HTCoefficients<T, Basis>& tree, const htucker::DimensionIndex& idx,
    const SepCoefficients<S, T, Index>& coeff);


template <typename T, SortingCriterion S, typename Index, typename Basis>
void
set_inplace(HTCoefficients<T, Basis>& tree, const htucker::DimensionIndex& idx,
            const SepCoefficients<S, T, Index>& coeff);


template <typename T, SortingCriterion S, typename Index, typename Basis>
void
axpy(HTCoefficients<T, Basis>& tree, const htucker::DimensionIndex& idx,
     const unsigned FLENS_DEFAULT_INDEXTYPE col, const T alpha,
     const Coefficients<S, T, Index>& coeff);


template <typename T, SortingCriterion S, typename Index, typename Basis>
void
axpy(HTCoefficients<T, Basis>& tree, const htucker::DimensionIndex& idx,
     const T alpha, const SepCoefficients<S, T, Index>& coeff);


template <typename T, SortingCriterion S, typename Index, typename Basis>
void
xpay(HTCoefficients<T, Basis>& tree, const htucker::DimensionIndex& idx,
     const unsigned FLENS_DEFAULT_INDEXTYPE col, const T alpha,
     const Coefficients<S, T, Index>& coeff);


template <typename T, SortingCriterion S, typename Index, typename Basis>
void
xpay(HTCoefficients<T, Basis>& tree, const htucker::DimensionIndex& idx,
     const T alpha, const SepCoefficients<S, T, Index>& coeff);


template <typename T, typename Basis>
Coefficients<Lexicographical, T, Index1D>
extract(const HTCoefficients<T, Basis>& tree,
        const htucker::DimensionIndex& idx,
        const unsigned FLENS_DEFAULT_INDEXTYPE col);


template <typename T, typename Basis>
SepCoefficients<Lexicographical, T, Index1D>
extract(const HTCoefficients<T, Basis>& tree,
        const htucker::DimensionIndex& idx);


template <typename T, typename Basis>
SepCoefficients<Lexicographical, T, Index1D>
extract(      HTCoefficients<T, Basis>& tree,
        const IndexSet<Index1D>& Lambda,
        const htucker::DimensionIndex& idx);


template <typename T, SortingCriterion S, typename Index>
SepCoefficients<S, T, Index>
add(const SepCoefficients<S, T, Index>& left,
    const SepCoefficients<S, T, Index>& right);


template <typename T, SortingCriterion S, typename Index>
SepCoefficients<S, T, Index>
operator+(const SepCoefficients<S, T, Index>& left,
          const SepCoefficients<S, T, Index>& right);


template <typename T, typename Optype>
SepCoefficients<Lexicographical, T, Index1D>
evalstandard(Sepop<Optype>& A,
             const SepCoefficients<Lexicographical, T, Index1D>& u,
             const std::size_t hashtablelength=193);


template <typename T, typename Optype>
SepCoefficients<Lexicographical, T, Index1D>
evalsimple(Sepop<Optype>& A,
           const SepCoefficients<Lexicographical, T, Index1D>& u,
           const std::size_t hashtablelength=193);


template <typename T, typename Optype>
SepCoefficients<Lexicographical, T, Index1D>
evallaplace(Sepop<Optype>& A,
            const SepCoefficients<Lexicographical, T, Index1D>& u,
            const std::size_t hashtablelength=193);


template <typename T, typename Optype>
SepCoefficients<Lexicographical, T, Index1D>
eval(Sepop<Optype>& A,
     const SepCoefficients<Lexicographical, T, Index1D>& u,
     const std::size_t hashtablelength=193);


template <typename T, typename Optype>
SepCoefficients<Lexicographical, T, Index1D>
operator*(Sepop<Optype>& A,
          const SepCoefficients<Lexicographical, T, Index1D>& u);


template <typename T, typename Optype>
SepCoefficients<Lexicographical, T, Index1D>
evalstandard(Sepop<Optype>& A,
             const SepCoefficients<Lexicographical, T, Index1D>& u,
             const std::vector<IndexSet<Index1D> >&  rows,
             const std::vector<IndexSet<Index1D> >&  cols,
             const std::size_t hashtablelength=193);


template <typename T, typename Optype>
SepCoefficients<Lexicographical, T, Index1D>
evalsimple(Sepop<Optype>& A,
           const SepCoefficients<Lexicographical, T, Index1D>& u,
           const std::vector<IndexSet<Index1D> >&  rows,
           const std::vector<IndexSet<Index1D> >&  cols,
           const std::size_t hashtablelength=193);


template <typename T, typename Optype>
SepCoefficients<Lexicographical, T, Index1D>
evallaplace(Sepop<Optype>& A,
            const SepCoefficients<Lexicographical, T, Index1D>& u,
            const std::vector<IndexSet<Index1D> >&  rows,
            const std::vector<IndexSet<Index1D> >&  cols,
            const std::size_t hashtablelength=193);


template <typename T, typename Optype>
SepCoefficients<Lexicographical, T, Index1D>
eval(Sepop<Optype>& A,
     const SepCoefficients<Lexicographical, T, Index1D>& u,
     const std::vector<IndexSet<Index1D> >&  rows,
     const std::vector<IndexSet<Index1D> >&  cols,
     const std::size_t hashtablelength=193);


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
evalstandard(Sepop<Optype>& A,
             const HTCoefficients<T, Basis>& u,
             const double eps = 1e-08,
             const std::size_t hashtablelength=193);


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
evalsimple(Sepop<Optype>& A,
           const HTCoefficients<T, Basis>& u,
           const double eps = 1e-08,
           const std::size_t hashtablelength=193);


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
evallaplace(Sepop<Optype>& A,
            const HTCoefficients<T, Basis>& u,
            const double eps = 1e-08,
            const std::size_t hashtablelength=193);


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
eval(Sepop<Optype>& A,
     const HTCoefficients<T, Basis>& u,
     const double eps = 1e-08,
     const std::size_t hashtablelength=193);


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
operator*(Sepop<Optype>& A,
          const HTCoefficients<T, Basis>& u);


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
evalstandard(Sepop<Optype>& A,
             const HTCoefficients<T, Basis>& u,
             const std::vector<IndexSet<Index1D> >&  rows,
             const std::vector<IndexSet<Index1D> >&  cols,
             const double eps = 1e-08,
             const std::size_t hashtablelength=193);


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
evalsimple_old(Sepop<Optype>& A,
               const HTCoefficients<T, Basis>& u,
               const std::vector<IndexSet<Index1D> >&  rows,
               const std::vector<IndexSet<Index1D> >&  cols,
               const double eps = 1e-08,
               const std::size_t hashtablelength=193);


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
evallaplace_old(Sepop<Optype>& A,
                const HTCoefficients<T, Basis>& u,
                const std::vector<IndexSet<Index1D> >&  rows,
                const std::vector<IndexSet<Index1D> >&  cols,
                const double eps = 1e-08,
                const std::size_t hashtablelength=193);


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
evalsimple(      Sepop<Optype>&            A,
                 HTCoefficients<T, Basis>& u,
           const std::vector<IndexSet<Index1D> >&  rows,
           const std::vector<IndexSet<Index1D> >&  cols,
           const std::size_t hashtablelength=193);


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
evallaplace(      Sepop<Optype>&            A,
                  HTCoefficients<T, Basis>& u,
            const std::vector<IndexSet<Index1D> >&  rows,
            const std::vector<IndexSet<Index1D> >&  cols,
            const std::size_t hashtablelength=193);


template <typename T, typename Optype, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
evallaplace(      Sepop<Optype>&            A,
            const flens::GeMatrix<
                  flens::FullStorage<T,
                  cxxblas::ColMajor> >&     U,
                  HTCoefficients<T, Basis>& u,
            const unsigned                    j,
            const IndexSet<Index1D>&        rows,
            const IndexSet<Index1D>&        cols,
            const std::size_t               hashtablelength=193);


template <typename T, typename Optype, typename Prec, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
evallaplace(      Sepop<Optype>&            A,
                  Prec&                     P,
            const flens::GeMatrix<
                  flens::FullStorage<T,
                  cxxblas::ColMajor> >&     U,
                  HTCoefficients<T, Basis>& u,
            const unsigned                  j,
            const IndexSet<Index1D>&        rows,
            const IndexSet<Index1D>&        cols,
            const std::size_t               hashtablelength=193);


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
eval(      Sepop<Optype>&            A,
           HTCoefficients<T, Basis>& u,
     const std::vector<IndexSet<Index1D> >&  rows,
     const std::vector<IndexSet<Index1D> >&  cols,
     const double eps = 1e-08,
     const std::size_t hashtablelength=193);


template <typename T, typename Optype, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
eval(      Sepop<Optype>&              A,
     const flens::GeMatrix<
           flens::FullStorage<T,
           cxxblas::ColMajor> >&       U,
           HTCoefficients<T, Basis>&   u,
     const unsigned                    j,
     const IndexSet<Index1D>&          rows,
     const IndexSet<Index1D>&          cols,
     const double                      eps = 1e-08,
     const std::size_t                 hashtablelength=193);


template <typename T, typename Optype, typename Prec, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
eval(      Sepop<Optype>&              A,
           Prec&                       P,
     const flens::GeMatrix<
           flens::FullStorage<T,
           cxxblas::ColMajor> >&       U,
           HTCoefficients<T, Basis>&   u,
     const unsigned                    j,
     const IndexSet<Index1D>&          rows,
     const IndexSet<Index1D>&          cols,
     const double                      eps = 1e-08,
     const std::size_t                 hashtablelength=193);


template <typename T, typename Optype, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
redeval(      Sepop<Optype>&              A,
        const flens::GeMatrix<
              flens::FullStorage<T,
              cxxblas::ColMajor> >&       Uj,
        const flens::GeMatrix<
              flens::FullStorage<T,
              cxxblas::ColMajor> >&       Pj,
              HTCoefficients<T, Basis>&   u,
        const unsigned                    j,
        const IndexSet<Index1D>&          rows,
        const IndexSet<Index1D>&          cols,
        const double                      eps = 1e-08,
        const std::size_t                 hashtablelength=193);


template <typename T, typename Optype, typename Prec, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
redeval_laplace(      Sepop<Optype>&              A,
                      Prec&                       P,
                const flens::GeMatrix<
                      flens::FullStorage<T,
                      cxxblas::ColMajor> >&       Uj,
                const flens::GeMatrix<
                      flens::FullStorage<T,
                      cxxblas::ColMajor> >&       Pj,
                      HTCoefficients<T, Basis>&   u,
                const unsigned                    j,
                const IndexSet<Index1D>&          rows,
                const IndexSet<Index1D>&          cols,
                const double                      eps = 1e-08,
                const std::size_t                 hashtablelength = 193);


template <typename T, typename Optype, typename Prec, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
redeval(      Sepop<Optype>&              A,
              Prec&                       P,
        const flens::GeMatrix<
              flens::FullStorage<T,
              cxxblas::ColMajor> >&       Uj,
        const flens::GeMatrix<
              flens::FullStorage<T,
              cxxblas::ColMajor> >&       Pj,
              HTCoefficients<T, Basis>&   u,
        const unsigned                    j,
        const IndexSet<Index1D>&          rows,
        const IndexSet<Index1D>&          cols,
        const double                      eps = 1e-08,
        const std::size_t                 hashtablelength = 193);


template <typename T, typename Optype, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
redeval(      Sepop<Optype>&                    A,
              Sepdiagscal<Basis>&               S,
        const flens::GeMatrix<
              flens::FullStorage<T,
              cxxblas::ColMajor> >&             Uj,
        const flens::GeMatrix<
              flens::FullStorage<T,
              cxxblas::ColMajor> >&             Pj,
              HTCoefficients<T, Basis>&         u,
        const unsigned                          j,
        const std::vector<IndexSet<Index1D> >&  rows,
        const std::vector<IndexSet<Index1D> >&  cols,
        const double                            eps = 1e-08,
        const std::size_t                       hashtablelength=193);


template <typename T, typename Optype, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
redeval(      Sepop<Optype>&              A,
              Sepdiagscal<Basis>&         S,
        const flens::GeMatrix<
              flens::FullStorage<T,
              cxxblas::ColMajor> >&       Uj,
        const flens::GeMatrix<
              flens::FullStorage<T,
              cxxblas::ColMajor> >&       Pj,
              HTCoefficients<T, Basis>&   u,
        const unsigned                    j,
        const IndexSet<Index1D>&          rows,
        const IndexSet<Index1D>&          cols,
        const double                      eps = 1e-08,
        const std::size_t                 hashtablelength = 193);


FLENS_DEFAULT_INDEXTYPE
maxlevel(const std::vector<IndexSet<Index1D> >& Lambda);


template <typename Basis>
typename Sepdiagscal<Basis>::T
compOmegamin2(const Basis& basis,
              const typename Sepdiagscal<Basis>::size_type d,
              const typename Sepdiagscal<Basis>::T order);


template <typename T>
T
compOmegamax2(const std::vector<IndexSet<Index1D> >& Lambda,
              const T order);


template <typename T>
T
compOmegamax4(const std::vector<IndexSet<Index1D> >& Lambda,
              const T order);


template <typename Basis>
typename Basis::T
compUnDistFac(const Basis& basis,
              const std::vector<IndexSet<Index1D> >& Lambda,
              const typename Basis::T order);


template <typename Basis>
typename Sepdiagscal<Basis>::T
compOmegamin4(const Basis& basis,
              const typename Sepdiagscal<Basis>::size_type d,
              const typename Sepdiagscal<Basis>::T order);


template <typename T, typename Basis>
T
compIndexscale(const SepCoefficients<Lexicographical, T, Index1D>& u,
               const Basis& basis, const T order);


template <typename T, typename Basis>
T
compIndexscale(const HTCoefficients<T, Basis>& u, const T order);


template <typename T, typename Basis>
T
compIndexscale2(const HTCoefficients<T, Basis>& u, const T order);


template <typename Basis>
typename Basis::T
compIndexscale2(const Basis& basis,
                const std::vector<IndexSet<Index1D> >& Lambda,
                const typename Basis::T order);


template <typename Basis>
typename Basis::T
compIndexscale(const Basis& basis,
               const std::vector<IndexSet<Index1D> >& Lambda,
               const typename Basis::T order);


template <typename Basis>
void
setScaling(Sepdiagscal<Basis>& S,
           const typename Sepdiagscal<Basis>::T eps);


template <typename T, typename Basis>
SepCoefficients<Lexicographical, T, Index1D>
eval(Sepdiagscal<Basis>& S,
     const SepCoefficients<Lexicographical, T, Index1D>& u);


template <typename T, typename Basis>
SepCoefficients<Lexicographical, T, Index1D>
operator*(Sepdiagscal<Basis>& S,
          const SepCoefficients<Lexicographical, T, Index1D>& u);


template <typename T, typename Basis>
SepCoefficients<Lexicographical, T, Index1D>
eval(Sepdiagscal<Basis>& S,
     const SepCoefficients<Lexicographical, T, Index1D>& u,
     const std::vector<IndexSet<Index1D> >&  cols);


template <typename T, typename Basis>
HTCoefficients<T, Basis>
eval(Sepdiagscal<Basis>& S,
     const HTCoefficients<T, Basis>& u,
     const double eps = 1e-08);


template <typename T, typename Basis>
HTCoefficients<T, Basis>
operator*(Sepdiagscal<Basis>& S,
          const HTCoefficients<T, Basis>& u);


template <typename T, typename Basis>
HTCoefficients<T, Basis>
eval(Sepdiagscal<Basis>& S,
     HTCoefficients<T, Basis>& u,
     const std::vector<IndexSet<Index1D> >&  cols,
     const double eps = 1e-08);


template <typename T, typename Basis>
HTCoefficients<T, Basis>
eval_notrunc(Sepdiagscal<Basis>&                     S,
             HTCoefficients<T, Basis>&               u,
             const std::vector<IndexSet<Index1D> >&  cols);


template <typename T, typename Basis>
HTCoefficients<T, Basis>
fixeval_notrunc(Sepdiagscal<Basis>&                      S,
                HTCoefficients<T, Basis>&                u,
                const std::vector<IndexSet<Index1D> >&   cols);


template <typename T, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
eval(      Sepdiagscal<Basis>&              S,
     const flens::GeMatrix
           <flens::FullStorage
           <T, cxxblas::ColMajor>>&         Uj,
           HTCoefficients<T, Basis>&        u,
     const unsigned                         j,
     const std::vector<IndexSet<Index1D> >& cols);


template <typename T, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
fixeval(      Sepdiagscal<Basis>&           S,
        const flens::GeMatrix
              <flens::FullStorage
              <T, cxxblas::ColMajor>>&      Uj,
              HTCoefficients<T, Basis>&     u,
     const unsigned                         j,
     const IndexSet<Index1D>&               cols);


template <typename Precon, typename T, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
prec(Precon&                                                     P,
     const
     flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& U,
     HTCoefficients<T, Basis>&                                   u,
     const unsigned                                              j,
     const IndexSet<Index1D>&                                    active);



template <typename Precon, typename T, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
prec(Precon&                                                     P,
     const flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& Pj,
     const
     flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& U,
     HTCoefficients<T, Basis>&                                   u,
     const unsigned                                              j,
     const IndexSet<Index1D>&                                    active);


template <typename Precon, typename T, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
precsq(Precon&                                                     P,
     const flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& Pj,
     const
     flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& U,
     HTCoefficients<T, Basis>&                                   u,
     const unsigned                                              j,
     const IndexSet<Index1D>&                                    active);


template <typename Precon, typename T, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
remove_prec(Precon&                                                     P,
            const
            flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& U,
            HTCoefficients<T, Basis>&                                   u,
            const unsigned                                              j,
            const IndexSet<Index1D>&                                    active);


template <typename Precon, typename T, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
precsq(Precon&                                                     P,
       const
       flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& U,
       HTCoefficients<T, Basis>&                                   u,
       const unsigned                                              j,
       const IndexSet<Index1D>&                                    active);


template <typename Precon, typename T, typename Basis>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
remove_precsq(Precon&                                      P,
              const
              flens::GeMatrix
              <flens::FullStorage<T, cxxblas::ColMajor> >& U,
              HTCoefficients<T, Basis>&                    u,
              const unsigned                               j,
              const IndexSet<Index1D>&                     active);


template <typename Precon, typename T, typename Basis>
void
rank1prec(Precon&                               P,
          HTCoefficients<T, Basis>&             u,
          const std::vector<IndexSet<Index1D>>& Lambda);


template <typename Precon, typename T, typename Basis>
void
rank1prec(Precon&                   P,
          HTCoefficients<T, Basis>& u,
          const unsigned            j,
          const IndexSet<Index1D>&  Lambda);


template <typename Precon, typename T, typename Basis>
void
remove_rank1prec(Precon&                               P,
                 HTCoefficients<T, Basis>&             u,
                 const std::vector<IndexSet<Index1D>>& Lambda);


template <typename Precon, typename T, typename Basis>
void
remove_rank1prec(Precon&                    P,
                 HTCoefficients<T, Basis>&  u,
                 const unsigned             j,
                 const IndexSet<Index1D>&   Lambda);


template <typename T, typename Basis>
HTCoefficients<T, Basis>
evalS2(Sepdiagscal<Basis>& S,
       HTCoefficients<T, Basis>& u,
       const std::vector<IndexSet<Index1D> >& cols,
       const double eps);


template <typename T, typename Basis>
void
assemble(Sepdiagscal<Basis>& S,
         HTCoefficients<T, Basis>& Stree,
         const std::vector<IndexSet<Index1D> >& cols);


template <typename Basis>
std::ostream& operator<<(std::ostream& s,
                         const Sepdiagscal<Basis>& S);


template <typename T, typename Basis, typename Index>
Coefficients<Lexicographical, T, Index>
contraction(      HTCoefficients<T, Basis>& u,
            const IndexSet<Index>& activex,
            const flens::DenseVector<flens::Array<T> >& sigmas,
            const FLENS_DEFAULT_INDEXTYPE dim);


template <typename T, typename Basis, typename Index>
void
contraction(      HTCoefficients<T, Basis>& u,
            const std::vector<IndexSet<Index> >& activex,
            const std::vector<flens::DenseVector<flens::Array<T> > >& sigmas,
            std::vector<Coefficients<Lexicographical, T, Index> >& ret);


template <typename T, typename Basis>
void
scal(const T alpha, HTCoefficients<T, Basis>& x);


template <typename T, typename Basis>
T
dot(const HTCoefficients<T, Basis>& x, const HTCoefficients<T, Basis>& y);


template <typename T, typename Basis>
T
nrm2(HTCoefficients<T, Basis>& x, bool isorth = false);


template <typename T, typename Basis>
void
axpy(const T alpha,
     const HTCoefficients<T, Basis>& x,
           HTCoefficients<T, Basis>& y);


template <typename T, typename Basis>
void
restrict(HTCoefficients<T, Basis>& f,
         const IndexSet<Index1D>& activex,
         const unsigned j);


template <typename T, typename Basis>
void
restrict(HTCoefficients<T, Basis>& f,
         const std::vector<IndexSet<Index1D> >& activex);


template <typename T, typename Basis>
void
extend(HTCoefficients<T, Basis>& f,
       const IndexSet<Index1D>& activex,
       const unsigned j);


template <typename T, typename Basis>
void
extend(HTCoefficients<T, Basis>& f,
       const std::vector<IndexSet<Index1D> >& activex);


template <typename T, typename Optype, typename Basis>
std::vector<flens::GeMatrix<
flens::FullStorage<T, cxxblas::ColMajor> >>
reduce_laplace(      Sepop<Optype>&                    A,
                     HTCoefficients<T, Basis>&         U,
               const std::vector<IndexSet<Index1D> >&  rows,
               const std::vector<IndexSet<Index1D> >&  cols,
               const std::size_t                       hashtablelength = 193);


template <typename T, typename Optype, typename Basis>
std::vector<flens::GeMatrix<
flens::FullStorage<T, cxxblas::ColMajor> >>
reduce(      Sepop<Optype>&                    A,
             HTCoefficients<T, Basis>&         U,
       const std::vector<IndexSet<Index1D> >&  rows,
       const std::vector<IndexSet<Index1D> >&  cols,
       const std::size_t                       hashtablelength = 193);


template <typename T, typename Basis>
std::vector<flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > >
reduce_rhs(const HTCoefficients<T, Basis>&         U,
           const HTCoefficients<T, Basis>&         b);

} // namespace lawa

#include <lawa/methods/adaptive/algorithms/coeffops.tcc>

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_COEFFOPS_H 1
