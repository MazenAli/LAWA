#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_SEPCOEFFICIENTS_TCC
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_SEPCOEFFICIENTS_TCC 1

#include <cassert>

namespace lawa
{

template <SortingCriterion S, typename T, typename Index>
SepCoefficients<S, T, Index>::SepCoefficients():
    rank_(0),
    dim_(0),
    coeffs(){}


template <SortingCriterion S, typename T, typename Index>
SepCoefficients<S, T, Index>::SepCoefficients(const size_type _rank,
                                              const size_type _dim):
    rank_(_rank),
    dim_(_dim),
    coeffs(rank()*dim()){}


template <SortingCriterion S, typename T, typename Index>
SepCoefficients<S, T, Index>::SepCoefficients(const CoeffVec& _coeffs,
                                              const size_type _rank,
                                              const size_type _dim):
    rank_(_rank),
    dim_(_dim),
    coeffs(_coeffs)
{
    assert(coeffs.size()==rank()*dim());
}


template <SortingCriterion S, typename T, typename Index>
void
SepCoefficients<S, T, Index>::resize(const size_type _rank,
                                     const size_type _dim)
{
    rank_ = _rank;
    dim_  = _dim;
    coeffs.resize(rank_*dim_);
}


template <SortingCriterion S, typename T, typename Index>
typename SepCoefficients<S, T, Index>::size_type
SepCoefficients<S, T, Index>::rank() const
{
    return rank_;
}


template <SortingCriterion S, typename T, typename Index>
typename SepCoefficients<S, T, Index>::size_type
SepCoefficients<S, T, Index>::dim() const
{
    return dim_;
}


template <SortingCriterion S, typename T, typename Index>
const typename SepCoefficients<S, T,Index>::CoeffVec&
SepCoefficients<S, T, Index>::getCoefficients() const
{
    return coeffs;
}


template <SortingCriterion S, typename T, typename Index>
typename SepCoefficients<S, T, Index>::CoeffVec&
SepCoefficients<S, T, Index>::getCoefficients()
{
    return coeffs;
}


template <SortingCriterion S, typename T, typename Index>
const typename SepCoefficients<S, T, Index>::Coeff&
SepCoefficients<S, T, Index>::getCoefficients(const size_type i,
                                              const size_type j) const
{
    assert(i>=1 && i<=rank());
    assert(j>=1 && j<=dim());
    return getCoefficients()[(j-1)*rank()+(i-1)];
}


template <SortingCriterion S, typename T, typename Index>
typename SepCoefficients<S, T, Index>::Coeff&
SepCoefficients<S, T, Index>::getCoefficients(const size_type i,
                                              const size_type j)
{
    assert(i>=1 && i<=rank());
    assert(j>=1 && j<=dim());
    return getCoefficients()[(j-1)*rank()+(i-1)];
}


template <SortingCriterion S, typename T, typename Index>
T
SepCoefficients<S, T, Index>::eval(const size_type i, const size_type j,
                                   const Index1D& index) const
{
    assert(i>=1 && i<=rank());
    assert(j>=1 && j<=dim());

    const auto end  = getCoefficients(i, j).end();
    const auto find = getCoefficients(i, j).find(index);
    if (find!=end) return getCoefficients(i, j).at(index);
    return (T) 0;
}


template <SortingCriterion S, typename T, typename Index>
T
SepCoefficients<S, T, Index>::eval(const IndexD& index) const
{
    assert(index.dim()==dim());
    T sum = 0.;
    for (size_type i=1; i<=rank(); ++i) {
        T prod = 1.;
        for (size_type j=1; j<=dim(); ++j) {
            prod *= eval(i, j, index(j));
        }
        sum += prod;
    }

    return sum;
}


template <SortingCriterion S, typename T, typename Index>
const typename SepCoefficients<S, T, Index>::Coeff&
SepCoefficients<S, T, Index>::operator()(const size_type i,
                                         const size_type j) const
{
    assert(i>=1 && i<=rank());
    assert(j>=1 && j<=dim());
    return getCoefficients(i, j);
}


template <SortingCriterion S, typename T, typename Index>
typename SepCoefficients<S, T, Index>::Coeff&
SepCoefficients<S, T, Index>::operator()(const size_type i,
                                         const size_type j)
{
    assert(i>=1 && i<=rank());
    assert(j>=1 && j<=dim());
    return getCoefficients(i, j);
}


template <SortingCriterion S, typename T, typename Index>
T
SepCoefficients<S, T, Index>::operator()(const size_type i,
                                         const size_type j,
                                         const Index1D& index) const
{
    assert(i>=1 && i<=rank());
    assert(j>=1 && j<=dim());
    return eval(i, j, index);
}


template <SortingCriterion S, typename T, typename Index>
T
SepCoefficients<S, T, Index>::operator()(const IndexD& index) const
{
    assert(index.dim()==dim());
    return eval(index);
}


template <SortingCriterion S, typename T, typename Index>
SepCoefficients<S, T, Index>&
SepCoefficients<S, T, Index>::
operator=(const SepCoefficients<S, T, Index>& copy)
{
    if ((rank()!=copy.rank()) || (dim()!=copy.dim())) {
        resize(copy.rank(), copy.dim());
    }

    for (size_type i=1; i<=copy.rank(); ++i) {
        for (size_type j=1; j<=copy.dim(); ++j) {
            getCoefficients(i, j) = copy(i, j);
        }
    }

    return *this;
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_SEPCOEFFICIENTS_TCC
