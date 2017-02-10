#include <cstddef>
#include <cassert>
#include <vector>
#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>


#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_COEFFFRAME_TCC
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_COEFFFRAME_TCC 1

namespace lawa
{


template <SortingCriterion S, typename T, typename _Index>
CoeffFrame<S, T, _Index>::CoeffFrame(const std::size_t _numCols_)
                                    :numCols_(_numCols_),
                                    U(numCols_)
{
    assert(numCols_ > 0);
}


template <SortingCriterion S, typename T, typename _Index>
CoeffFrame<S, T, _Index>::CoeffFrame(const CoeffFrame<S, T, _Index>& copy)
                                    :numCols_(copy.numCols()),
                                    U(copy.getFrame())
{
    assert(numCols_ > 0);
}


template <SortingCriterion S, typename T, typename _Index>
const std::vector<Coefficients<S, T, _Index> >&
CoeffFrame<S, T, _Index>::getFrame() const
{
    return U;
}


template <SortingCriterion S, typename T, typename _Index>
std::size_t
CoeffFrame<S, T, _Index>::numCols() const
{
    return numCols_;
}


template <SortingCriterion S, typename T, typename _Index>
const Coefficients<S, T, _Index>&
CoeffFrame<S, T, _Index>::operator[] (const std::size_t col_num) const
{
    assert(col_num >= 1 && col_num <= numCols_);
    return U[col_num-1];
}


template <SortingCriterion S, typename T, typename _Index>
T
CoeffFrame<S, T, _Index>::operator() (const _Index& lambda, const std::size_t col_num) const
{
    assert(col_num >= 1 && col_num <= numCols_);
    return U[col_num-1][lambda];
}


template <SortingCriterion S, typename T, typename _Index>
Coefficients<S, T, _Index>&
CoeffFrame<S, T, _Index>::operator[] (const std::size_t col_num)
{
    assert(col_num >= 1 && col_num <= numCols_);
    return U[col_num-1];
}


template <SortingCriterion S, typename T, typename _Index>
T&
CoeffFrame<S, T, _Index>::operator() (const _Index& lambda, const std::size_t col_num)
{
    assert(col_num >= 1 && col_num <= numCols_);
    return U[col_num-1][lambda];
}


} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_COEFFFRAME_TCC
