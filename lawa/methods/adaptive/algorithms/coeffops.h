#ifndef LAWA_METHODS_ADAPTIVE_ALGORITHMS_COEFFOPS_H
#define LAWA_METHODS_ADAPTIVE_ALGORITHMS_COEFFOPS_H 1

#include <iostream>
#include <vector>
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


template <SortingCriterion S, typename T, typename Index>
std::ostream& operator<<(std::ostream& s,
                         const SepCoefficients<S, T, Index>& coeffs);


template <SortingCriterion S, typename T, typename Index, typename Basis>
unsigned long
maxintind(const Coefficients<S, T, Index>& coeffs, const Basis& basis);


template <typename T, SortingCriterion S, typename Index, typename Basis>
void
set(HTCoefficients<T, Basis>& tree,
    const SepCoefficients<S, T, Index>& cp);


template <typename T, SortingCriterion S, typename Index, typename Basis>
void
set(HTCoefficients<T, Basis>& tree, const htucker::DimensionIndex& idx,
    const unsigned long col, const Coefficients<S, T, Index>& coeff);


template <typename T, SortingCriterion S, typename Index, typename Basis>
void
set(HTCoefficients<T, Basis>& tree, const htucker::DimensionIndex& idx,
    const SepCoefficients<S, T, Index>& coeff);


template <typename T, SortingCriterion S, typename Index, typename Basis>
void
axpy(HTCoefficients<T, Basis>& tree, const htucker::DimensionIndex& idx,
     const unsigned long col, const T alpha,
     const Coefficients<S, T, Index>& coeff);


template <typename T, SortingCriterion S, typename Index, typename Basis>
void
axpy(HTCoefficients<T, Basis>& tree, const htucker::DimensionIndex& idx,
     const T alpha, const SepCoefficients<S, T, Index>& coeff);


template <typename T, SortingCriterion S, typename Index, typename Basis>
void
xpay(HTCoefficients<T, Basis>& tree, const htucker::DimensionIndex& idx,
     const unsigned long col, const T alpha,
     const Coefficients<S, T, Index>& coeff);


template <typename T, SortingCriterion S, typename Index, typename Basis>
void
xpay(HTCoefficients<T, Basis>& tree, const htucker::DimensionIndex& idx,
     const T alpha, const SepCoefficients<S, T, Index>& coeff);


template <typename T, typename Basis>
Coefficients<Lexicographical, T, Index1D>
extract(const HTCoefficients<T, Basis>& tree,
        const htucker::DimensionIndex& idx,
        const unsigned long col);


template <typename T, typename Basis>
SepCoefficients<Lexicographical, T, Index1D>
extract(const HTCoefficients<T, Basis>& tree,
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
evalsimple(Sepop<Optype>& A,
           const HTCoefficients<T, Basis>& u,
           const std::vector<IndexSet<Index1D> >&  rows,
           const std::vector<IndexSet<Index1D> >&  cols,
           const std::size_t hashtablelength=193);


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
evallaplace(Sepop<Optype>& A,
            const HTCoefficients<T, Basis>& u,
            const std::vector<IndexSet<Index1D> >&  rows,
            const std::vector<IndexSet<Index1D> >&  cols,
            const std::size_t hashtablelength=193);


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
eval(Sepop<Optype>& A,
     const HTCoefficients<T, Basis>& u,
     const std::vector<IndexSet<Index1D> >&  rows,
     const std::vector<IndexSet<Index1D> >&  cols,
     const double eps = 1e-08,
     const std::size_t hashtablelength=193);


template <typename T, typename Basis>
int
maxlevel(const HTCoefficients<T, Basis>& u);


template <typename Basis>
typename Sepdiagscal<Basis>::T
compOmegamin2(const Basis& basis,
              const typename Sepdiagscal<Basis>::size_type d,
              const typename Sepdiagscal<Basis>::T order);


template <typename T, typename Basis>
T
compOmegamax2(const HTCoefficients<T, Basis>& u, const T order);


template <typename T, typename Basis>
T
compUnDistFac(const HTCoefficients<T, Basis>& u, const T order);


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
     const HTCoefficients<T, Basis>& u,
     const std::vector<IndexSet<Index1D> >&  cols,
     const double eps = 1e-08);


template <typename T, typename Basis>
HTCoefficients<T, Basis>
evalS2(Sepdiagscal<Basis>& S,
       const HTCoefficients<T, Basis>& u,
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
contraction(const HTCoefficients<T, Basis>& u,
            const IndexSet<Index>& activex,
            const flens::DenseVector<flens::Array<T> >& sigmas,
            const int dim);


template <typename T, typename Basis, typename Index>
void
contraction(const HTCoefficients<T, Basis>& u,
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

} // namespace lawa

#include <lawa/methods/adaptive/algorithms/coeffops.tcc>

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_COEFFOPS_H 1
