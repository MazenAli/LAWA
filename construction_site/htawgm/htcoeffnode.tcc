#include <iostream>
#include <cstddef>
#include <cassert>
#include <stdlib.h>
#include <flens/flens.cxx>
#include <htucker/htucker.h>
#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <construction_site/coeffframe.h>

#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_HTCOEFFNODE_TCC
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_HTCOEFFNODE_TCC 1

namespace lawa
{


template <typename T, typename _Basis, typename _Index>
HTCoeffNode<T, _Basis, _Index>::HTCoeffNode
                    (htucker::HTuckerTreeNode<T>& _htnode,
                    const _Basis& _basis)
                    :htnode(_htnode),
                    basis(_basis),
                    numCols_(htnode.getUorB().numCols()){}


template <typename T, typename _Basis, typename _Index>
HTCoeffNode<T, _Basis, _Index>::HTCoeffNode
                    (const HTCoeffNode<T, _Basis, _Index>& copy)
                    :htnode(const_cast  <htucker::HTuckerTreeNode<T>&>
                                        (copy.getNode())),
                    basis(copy.getBasis()),
                    numCols_(0)
{
    std::cerr <<
"error: HTCoeffNode copy constructor being called!"
    << std::endl;
    exit(EXIT_FAILURE);
}


template <typename T, typename _Basis, typename _Index>
const htucker::HTuckerTreeNode<T>&
HTCoeffNode<T, _Basis, _Index>::getNode() const
{
    return htnode;
}


template <typename T, typename _Basis, typename _Index>
const IndexSet<_Index>&
HTCoeffNode<T, _Basis, _Index>::getActivex(const std::size_t col_num) const
{
    assert(activex.size()>0);
    assert(col_num >= 1 && col_num <= numCols_);
    return activex[col_num-1];
}


template <typename T, typename _Basis, typename _Index>
std::size_t
HTCoeffNode<T, _Basis, _Index>::numCols() const
{
    return numCols_;
}


template <typename T, typename _Basis, typename _Index>
const _Basis&
HTCoeffNode<T, _Basis, _Index>::getBasis() const
{
    return basis;
}


template <typename T, typename _Basis, typename _Index>
template <SortingCriterion S>
void
HTCoeffNode<T, _Basis, _Index>::addCoeff(const CoeffFrame<S, T, _Index>& u)
{
    if (activex.size() == 0) {
        std::cerr   << "HTCoeffNode must be initiliazed via setCoeff!"
                    << std::endl;
        exit(EXIT_FAILURE);
    }

    using flens::_;
    typedef typename CoeffFrame<S, T, _Index>::const_iterator   const_it;
    typedef flens::GeMatrix<flens::FullStorage<T,
                            flens::ColMajor>>                   GeMatrix;

    unsigned long numRows = htnode.getUorB().numRows();
    unsigned long current(0), max(0);

    assert(numCols_ > 0);
    assert(numCols_ == u.numCols());
    for (std::size_t i = 1; i <= numCols_; ++i) {
        for (const_it it = u[i].cbegin(); it != u[i].cend(); ++it) {
            activex[i-1].insert((*it).first);
            current = mapCoeff((*it).first, basis);
            if (current > max) {
                max = current;
            }
        }
    }

    GeMatrix& U = const_cast<GeMatrix&>(htnode.getUorB()); // in-place
    if (max > numRows && numRows > 0) {
        GeMatrix copy(U);
        U.resize(max, numCols_);
        U(_(1,numRows), _) = copy;
    } else if (numRows == 0 && max > 0) {
        U.resize(max, numCols_);
    } else if (numRows == 0 && max == 0) {
        U.resize(1, numCols_);
        return;
    }

    for (std::size_t i = 1; i <= numCols_; ++i) {
        for (const_it it = u[i].cbegin(); it != u[i].cend(); ++it) {
            current = mapCoeff((*it).first, basis);
            U(current, i) += (*it).second;
        }
    }
}


template <typename T, typename _Basis, typename _Index>
template <SortingCriterion S>
void
HTCoeffNode<T, _Basis, _Index>::setCoeff(const CoeffFrame<S, T, _Index>& u)
{
    typedef typename CoeffFrame<S, T, _Index>::const_iterator   const_it;
    typedef flens::GeMatrix<flens::FullStorage<T,
                            cxxblas::ColMajor>>                 GeMatrix;
    numCols_ = u.numCols();
    assert(numCols_ > 0);
    activex.clear();
    activex.resize(numCols_);

    unsigned long numRows(0), current(0);
    for (std::size_t i = 1; i <= numCols_; ++i) {
        for (const_it it = u[i].cbegin(); it != u[i].cend(); ++it) {
            activex[i-1].insert((*it).first);
            current = mapCoeff((*it).first, basis);
            if (current > numRows) {
                numRows = current;
            }
        }
    }

    GeMatrix& U = const_cast<GeMatrix&>(htnode.getUorB()); // in-place
    if (numRows == 0) {
        U.resize(1, numCols_);
        return;
    }

    U.resize(numRows, numCols_);
    for (std::size_t i = 1; i <= numCols_; ++i) {
        for (const_it it = u[i].cbegin(); it != u[i].cend(); ++it) {
            current = mapCoeff((*it).first, basis);
            U(current, i) = (*it).second;
        }
    }
}


template <typename T, typename _Basis, typename _Index>
const flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor>>&
HTCoeffNode<T, _Basis, _Index>::getFrame() const
{
    return htnode.getUorB();
}


template<typename T, typename _Basis, typename _Index>
void
HTCoeffNode<T, _Basis, _Index>::printActive()
{
    std::cout   << "-----   Printing active wavelet indices -----"
                << std::endl;
    for (int i = 0; i < numCols_; ++i)
    {
        std::cout << "Column number " << i+1 << std::endl;
        std::cout << activex[i] << std::endl;
    }
    std::cout   << "-----               Done                -----"
                << std::endl;
}


template <typename T, typename _Basis, typename _Index>
T
HTCoeffNode<T, _Basis, _Index>::operator() (const _Index& lambda,
                                            const std::size_t col_num) const
{
    assert(activex.size()>0);
    assert(col_num >= 1 && col_num <= numCols_);

    typename IndexSet<_Index>:: const_iterator
                                found = activex[col_num-1].find(lambda);
    if (found == activex[col_num-1].end()) {
        return 0;
    }

    return getFrame()(mapCoeff(lambda, basis), col_num);
}


template <typename T, typename _Basis, typename _Index>
bool
HTCoeffNode<T, _Basis, _Index>::isActive(   const _Index& lambda,
                                            const std::size_t col_num) const
{
    assert(activex.size()>0);
    assert(col_num >= 1 && col_num <= numCols_);
    typename IndexSet<_Index>:: const_iterator
                                found = activex[col_num-1].find(lambda);

    if (found == activex[col_num-1].end()) {
        return false;
    }
    return true;
}


template <typename T, typename _Basis, typename _Index>
typename HTCoeffNode<T, _Basis, _Index>::iterator
HTCoeffNode<T, _Basis, _Index>::begin(const std::size_t j)
{
    assert(activex.size() > 0);
    assert(j >= 1 && j <= numCols_);
    return activex[j-1].begin();
}


template <typename T, typename _Basis, typename _Index>
typename HTCoeffNode<T, _Basis, _Index>::const_iterator
HTCoeffNode<T, _Basis, _Index>::cbegin(const std::size_t j) const
{
    assert(activex.size() > 0);
    assert(j >= 1 && j <= numCols_);
    return activex[j-1].cbegin();
}


template <typename T, typename _Basis, typename _Index>
typename HTCoeffNode<T, _Basis, _Index>::iterator
HTCoeffNode<T, _Basis, _Index>::end(const std::size_t j)
{
    assert(activex.size()>0);
    assert(j >= 1 && j <= numCols_);
    return activex[j-1].end();
}


template <typename T, typename _Basis, typename _Index>
typename HTCoeffNode<T, _Basis, _Index>::const_iterator
HTCoeffNode<T, _Basis, _Index>::cend(const std::size_t j) const
{
    assert(activex.size()>0);
    assert(j >= 1 && j <= numCols_);
    return activex[j-1].cend();
}


template <typename _Index, typename _Basis>
unsigned long
mapCoeff(const _Index&, const _Basis&)
{
    std::cerr   << "error: mapCoeff not implemented for given Index type"
                << std::endl;
    exit(EXIT_FAILURE);
}


template <typename _Basis>
unsigned long
mapCoeff(const Index1D& lambda, const _Basis& basis)
{
    long j       = lambda.j;
    long k       = lambda.k;
    XType type   = lambda.xtype;
    long j0      = basis.j0;
    long offset  = 0;

    assert(type == XBSpline || type == XWavelet);
    if (type == XBSpline) {
        assert(j == j0);
        assert( k >= basis.mra.rangeI(j).firstIndex() &&
                k <= basis.mra.rangeI(j).lastIndex());
        offset -= basis.mra.rangeI(j).firstIndex();
        ++offset;
        return k + offset;
    } else {
        assert(j >= j0);
        assert( k >= basis.rangeJ(j).firstIndex() &&
                k <= basis.rangeJ(j).lastIndex());
        offset += basis.mra.cardI(j);
        offset -= basis.rangeJ(j).firstIndex();
        ++offset;
        return k + offset;
    }
}


} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_HTCOEFFNODE_TCC
