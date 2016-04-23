#ifndef LAWA_METHODS_ADAPTIVE_OPERATORS_OPERATORSD_SEPOP_TCC
#define LAWA_METHODS_ADAPTIVE_OPERATORS_OPERATORSD_SEPOP_TCC 1

#include <cassert>
#include <iostream>
#include <htucker/htucker.h>
#include <lawa/methods/adaptive/algorithms/coeffops.h>
#include <lawa/methods/adaptive/datastructures/treecoefficients1d.h>

namespace lawa
{

template <typename Optype>
Sepop<Optype>::Sepop(const Opvec& _ops,
                     const size_type _rank, const size_type _dim):
    ops_(_ops),
    rank_(_rank),
    dim_(_dim),
    indrows(_dim),
    indcols(_dim)
{
    assert(ops_.size()==rank()*dim() ||
           (ops_.size()==rank() && rank()<=dim()));

    if (ops_.size()==rank()*dim()) {
        type_ = standard;
    } else {
        type_ = simple;
    }
}


template <typename Optype>
Sepop<Optype>::Sepop(Optype& op,
                     const size_type _rank, const size_type _dim):
    ops_(1, &op),
    rank_(_rank),
    dim_(_dim),
    indrows(_dim),
    indcols(_dim),
    type_(laplace)
{
    assert(rank()<=dim());
}


template <typename Optype>
typename Sepop<Optype>::size_type
Sepop<Optype>::rank() const
{
    return rank_;
}


template <typename Optype>
typename Sepop<Optype>::size_type
Sepop<Optype>::dim() const
{
    return dim_;
}


template <typename Optype>
SepopType
Sepop<Optype>::type() const
{
    return type_;
}


template <typename Optype>
const typename Sepop<Optype>::IndexSetVec&
Sepop<Optype>::getrows() const
{
    return indrows;
}


template <typename Optype>
const typename Sepop<Optype>::IndexSetVec&
Sepop<Optype>::getcols() const
{
    return indcols;
}


template <typename Optype>
const IndexSet<Index1D>&
Sepop<Optype>::getrows(const size_type j) const
{
    assert(j>=1 && j<=dim());
    return indrows[j-1];
}


template <typename Optype>
const IndexSet<Index1D>&
Sepop<Optype>::getcols(const size_type j) const
{
    assert(j>=1 && j<=dim());
    return indcols[j-1];
}


template <typename Optype>
void
Sepop<Optype>::setrows(const IndexSetVec& _indrows)
{
    assert(dim()==_indrows.size());
    indrows = _indrows;
}


template <typename Optype>
void
Sepop<Optype>::setcols(const IndexSetVec& _indcols)
{
    assert(dim()==_indcols.size());
    indcols = _indcols;
}


template <typename Optype>
void
Sepop<Optype>::setrows(const IndexSet<Index1D>& _indrows, const size_type j)
{
    assert(j>=1 && j<=dim());
    indrows[j-1] = _indrows;
}


template <typename Optype>
void
Sepop<Optype>::setcols(const IndexSet<Index1D>& _indcols, const size_type j)
{
    assert(j>=1 && j<=dim());
    indcols[j-1] = _indcols;
}


template <typename Optype>
const typename Sepop<Optype>::Opvec&
Sepop<Optype>::ops() const
{
    return ops_;
}


template <typename Optype>
typename Sepop<Optype>::Opvec&
Sepop<Optype>::ops()
{
    return ops_;
}


template <typename Optype>
const Optype&
Sepop<Optype>::ops(const size_type i, const size_type j) const
{
    if (type()==standard) {
        assert(i>=1 && i<=rank());
        assert(j>=1 && j<=dim());
    } else if (type()==simple) {
        assert(i>=1 && i<=rank());
        assert(j==1);
    } else {
        assert(i==1 && j==1);
    }

    return *ops()[(j-1)*rank()+(i-1)];

}


template <typename Optype>
Optype&
Sepop<Optype>::ops(const size_type i, const size_type j)
{
    if (type()==standard) {
        assert(i>=1 && i<=rank());
        assert(j>=1 && j<=dim());
    } else if (type()==simple) {
        assert(i>=1 && i<=rank());
        assert(j==1);
    } else {
        assert(i==1 && j==1);
    }

    return *ops()[(j-1)*rank()+(i-1)];
}


template <typename Optype>
const Optype&
Sepop<Optype>::operator()(const size_type i, const size_type j) const
{
    if (type()==standard) {
        assert(i>=1 && i<=rank());
        assert(j>=1 && j<=dim());
    } else if (type()==simple) {
        assert(i>=1 && i<=rank());
        assert(j==1);
    } else {
        assert(i==1 && j==1);
    }

    return ops(i, j);
}


template <typename Optype>
Optype&
Sepop<Optype>::operator()(const size_type i, const size_type j)
{
    if (type()==standard) {
        assert(i>=1 && i<=rank());
        assert(j>=1 && j<=dim());
    } else if (type()==simple) {
        assert(i>=1 && i<=rank());
        assert(j==1);
    } else {
        assert(i==1 && j==1);
    }

    return ops(i, j);
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_OPERATORS_OPERATORSD_SEPOP_TCC
