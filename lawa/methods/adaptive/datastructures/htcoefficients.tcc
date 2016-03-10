#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_HTCOEFFICIENTS_TCC
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_HTCOEFFICIENTS_TCC 1

#include <lawa/methods/adaptive/datastructures/indexops.h>
#include <flens/flens.cxx>
#include <cassert>

namespace lawa
{

template <typename T, typename Basis>
HTCoefficients<T, Basis>::HTCoefficients(const int d, const Basis& _basis):
    httree(d),
    basis_(&_basis){}


template <typename T, typename Basis>
HTCoefficients<T, Basis>::HTCoefficients(const int d, const double split,
                                  const Basis& _basis):
    httree(d, split),
    basis_(&_basis){}


template <typename T, typename Basis>
int
HTCoefficients<T, Basis>::dim() const
{
    return httree.dim();
}


template <typename T, typename Basis>
void
HTCoefficients<T, Basis>::print() const
{
    httree.print();
}


template <typename T, typename Basis>
const typename HTCoefficients<T, Basis>::HTTree&
HTCoefficients<T, Basis>::tree() const
{
    return httree;
}


template <typename T, typename Basis>
typename HTCoefficients<T, Basis>::HTTree&
HTCoefficients<T, Basis>::tree()
{
    return httree;
}


template <typename T, typename Basis>
const Basis&
HTCoefficients<T, Basis>::basis() const
{
    return *basis_;
}


template <typename T, typename Basis>
void
HTCoefficients<T, Basis>::orthogonolize()
{
    tree().orthogonolize();
}


template <typename T, typename Basis>
void
HTCoefficients<T, Basis>::truncate(const int rank, bool isorth)
{
    tree().truncate(rank, isorth);
}


template <typename T, typename Basis>
void
HTCoefficients<T, Basis>::truncate(double eps, bool isorth)
{
    tree().truncate(eps, isorth);
}


template <typename T, typename Basis>
T
HTCoefficients<T, Basis>::eval(const IndexD& index) const
{
    assert(index.dim()==(unsigned) dim());

    typedef htucker::DimensionIndex                          IDX;
    typedef typename flens::DenseVector<flens::Array<int> >  IDV;

    IDX idx(dim());
    IDV intindex(dim());
    for (int i=1; i<=dim(); ++i) {
        intindex(i) = maptoint(index(i), basis());
    }
    idx.setValue(intindex);

    return httree.evaluate(idx);
}


template <typename T, typename Basis>
T
HTCoefficients<T, Basis>::eval(const IndexD& index, const int vardim) const
{
    assert(index.dim()==(unsigned) dim());

    htucker::DimensionIndex                 idx(dim());
    flens::DenseVector<flens::Array<int> >  intindex(dim());
    for (int i=1; i<=dim(); ++i) {
        intindex(i) = maptoint(index(i), basis());
    }
    idx.setValue(intindex);

    return httree.vec_evaluate(idx, vardim);
}


template <typename T, typename Basis>
T
HTCoefficients<T, Basis>::operator()(const IndexD& index) const
{
    assert(index.dim()==(unsigned) dim());
    return eval(index);
}


template <typename T, typename Basis>
T
HTCoefficients<T, Basis>::operator()(const IndexD& index,
                                     const int vardim) const
{
    assert(index.dim()==(unsigned) dim());
    return eval(index, vardim);
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_HTCOEFFICIENTS_TCC
