#ifndef LAWA_RIGHTHANDSIDES_SEPARABLERHSD_TCC
#define LAWA_RIGHTHANDSIDES_SEPARABLERHSD_TCC 1

#include <cassert>
#include <utility>

namespace lawa
{

template <typename T, typename Basis>
SeparableRHSD<T, Basis>::SeparableRHSD(const Basis&     _basis,
                                       const FuncType&  _F,
                                       const FLENS_DEFAULT_INDEXTYPE        order):
    rank_(_F.rank()),
    dim_(_F.dim()),
    basis(&_basis),
    F(&_F),
    deltas(nullptr),
    derivs(nullptr),
    data(rank_*dim_)
{
    for (size_type j=1; j<=dim(); ++j) {
        for (size_type i=1; i<=rank(); ++i) {
            integral.push_back(IntegralType((*F)(i, j), *basis));
            getInt(i, j).quadrature.setOrder(order);
        }
    }
}


template <typename T, typename Basis>
SeparableRHSD<T, Basis>::SeparableRHSD(const Basis&    _basis,
                                       const FuncType& _F,
                                       const Matvec&   _deltas,
                                       const IntMat&   _derivs,
                                       const FLENS_DEFAULT_INDEXTYPE       order):
    rank_(_F.rank()),
    dim_(_F.dim()),
    basis(&_basis),
    F(&_F),
    deltas(&_deltas),
    derivs(&_derivs),
    data(rank_*dim_)
{
    assert(!deltas->size() || deltas->size()==dim()*rank());
    assert((!derivs->numRows() && !derivs->numCols())
           || ((size_type) derivs->numRows()==rank() &&
               (size_type) derivs->numCols()==dim()));

    for (size_type j=1; j<=dim(); ++j) {
        for (size_type i=1; i<=rank(); ++i) {
            integral.push_back(IntegralType((*F)(i,j), *basis));
            getInt(i, j).quadrature.setOrder(order);
        }
    }
}


template <typename T, typename Basis>
typename SeparableRHSD<T, Basis>::size_type
SeparableRHSD<T, Basis>::rank() const
{
    return rank_;
}


template <typename T, typename Basis>
typename SeparableRHSD<T, Basis>::size_type
SeparableRHSD<T, Basis>::dim() const
{
    return dim_;
}


template <typename T, typename Basis>
const typename SeparableRHSD<T, Basis>::IntegralType&
SeparableRHSD<T, Basis>::getInt(const size_type i, const size_type j) const
{
    assert(i>=1 && i<=rank());
    assert(j>=1 && j<=dim());
    return integral[(j-1)*rank()+(i-1)];
}


template <typename T, typename Basis>
typename SeparableRHSD<T, Basis>::IntegralType&
SeparableRHSD<T, Basis>::getInt(const size_type i, const size_type j)
{
    assert(i>=1 && i<=rank());
    assert(j>=1 && j<=dim());
    return integral[(j-1)*rank()+(i-1)];
}


template <typename T, typename Basis>
const typename SeparableRHSD<T, Basis>::GeMat&
SeparableRHSD<T, Basis>::getDeltas(const size_type i, const size_type j) const
{
    assert(i>=1 && i<=rank());
    assert(j>=1 && j<=dim());
    assert(deltas);
    return (*deltas)[(j-1)*rank()+(i-1)];
}


template <typename T, typename Basis>
T
SeparableRHSD<T, Basis>::eval(const size_type i, const size_type j,
                              const Index1D& index)
{
    assert(i>=1 && i<=rank());
    assert(j>=1 && j<=dim());
    assert(data.size() == rank()*dim());

    // Typedefs
    typedef typename GeMat::IndexType IndexType;

    // Check if in data
    if (data[(j-1)*rank()+(i-1)].find(index) !=
        data[(j-1)*rank()+(i-1)].end())
        return data[(j-1)*rank()+(i-1)][index];

    // Compute otherwise
    int derivij = 0;
    if (derivs) {
        if (derivs->numRows() && derivs->numCols()) {
            derivij = (*derivs)(i, j);
        }
    }

    T ret = getInt(i, j)(index.j, index.k, index.xtype, derivij);
    if (deltas) {
        if (deltas->size()) {
            const GeMat& D = getDeltas(i, j);
            for (IndexType i=D.firstRow();
                 i<=D.lastRow(); ++i) {
                ret += D(i, 2)*(*basis).generator(index.xtype)(D(i, 1),
                                                            index.j, index.k,
                                                            derivij);
            }
        }
    }

    data[(j-1)*rank()+(i-1)][index] = ret;

    return ret;
}


template <typename T, typename Basis>
typename SeparableRHSD<T, Basis>::Coeff1D
SeparableRHSD<T, Basis>::eval(const size_type i, const size_type j,
                              const IndexSet<Index1D>& indexset)
{
    assert(i>=1 && i<=rank());
    assert(j>=1 && j<=dim());

    Coeff1D ret;
    T val;
    for (auto& lambda : indexset) {
        val = eval(i, j, lambda);
        ret[lambda] = val;
    }

    return ret;
}


template <typename T, typename Basis>
T
SeparableRHSD<T, Basis>::eval(const IndexD& index)
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


template <typename T, typename Basis>
T
SeparableRHSD<T, Basis>::operator()(const size_type i, const size_type j,
                                    const Index1D& index)
{
    assert(i>=1 && i<=rank());
    assert(j>=1 && j<=dim());
    return eval(i, j, index);
}


template <typename T, typename Basis>
typename SeparableRHSD<T, Basis>::Coeff1D
SeparableRHSD<T, Basis>::operator()(const size_type i, const size_type j,
                                    const IndexSet<Index1D>& indexset)
{
    assert(i>=1 && i<=rank());
    assert(j>=1 && j<=dim());
    return eval(i, j, indexset);
}


template <typename T, typename Basis>
T
SeparableRHSD<T, Basis>::operator()(const IndexD& index)
{
    assert(index.dim()==dim());
    return eval(index);
}

} // namespace lawa

#endif // LAWA_RIGHTHANDSIDES_SEPARABLERHSD_TCC
