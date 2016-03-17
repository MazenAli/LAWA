#ifndef LAWA_METHODS_ADAPTIVE_ALGORITHMS_COEFFOPS_H
#define LAWA_METHODS_ADAPTIVE_ALGORITHMS_COEFFOPS_H 1

#include <iostream>
#include <vector>
#include <lawa/settings/enum.h>
#include <lawa/righthandsides/separablerhsd.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/datastructures/sepcoefficients.h>
#include <lawa/methods/adaptive/datastructures/htcoefficients.h>
#include <lawa/methods/adaptive/operators/operatorsd/sepop.h>

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
eval(Sepop<Optype>& A,
     const SepCoefficients<Lexicographical, T, Index1D>& u,
     const std::size_t hashtablelength=193);


template <typename T, typename Optype>
SepCoefficients<Lexicographical, T, Index1D>
operator*(Sepop<Optype>& A,
          const SepCoefficients<Lexicographical, T, Index1D>& u);


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
eval(Sepop<Optype>& A,
     const HTCoefficients<T, Basis>& u,
     const std::size_t hashtablelength=193);


template <typename T, typename Optype, typename Basis>
HTCoefficients<T, Basis>
operator*(Sepop<Optype>& A,
          const HTCoefficients<T, Basis>& u);

} // namespace lawa

#include <lawa/methods/adaptive/algorithms/coeffops.tcc>

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_COEFFOPS_H 1
