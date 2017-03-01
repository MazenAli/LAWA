#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_HTCOEFFICIENTS_TCC
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_HTCOEFFICIENTS_TCC 1

#include <lawa/methods/adaptive/algorithms/indexops.h>
#include <flens/flens.cxx>
#include <cassert>
#include <cstddef>
#include <limits>

namespace lawa
{

template <typename T, typename Basis>
HTCoefficients<T, Basis>::HTCoefficients():
    httree(),
    basis_(nullptr),
    map_(nullptr){}


template <typename T, typename Basis>
HTCoefficients<T, Basis>::HTCoefficients(const FLENS_DEFAULT_INDEXTYPE d,
                                         const Basis& _basis,
                                         Maptype& _map):
    httree(d),
    basis_(&_basis),
    map_(&_map){}


template <typename T, typename Basis>
HTCoefficients<T, Basis>::HTCoefficients(const FLENS_DEFAULT_INDEXTYPE d,
                                         const double split,
                                         const Basis& _basis,
                                         Maptype& _map):
    httree(d, split),
    basis_(&_basis),
    map_(&_map){}


template <typename T, typename Basis>
FLENS_DEFAULT_INDEXTYPE
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
    assert(basis_);
    return *basis_;
}


template <typename T, typename Basis>
const typename HTCoefficients<T, Basis>::Maptype&
HTCoefficients<T, Basis>::map() const
{
    assert(map_);
    return *map_;
}


template <typename T, typename Basis>
typename HTCoefficients<T, Basis>::Maptype&
HTCoefficients<T, Basis>::map()
{
    assert(map_);
    return *map_;
}


template <typename T, typename Basis>
void
HTCoefficients<T, Basis>::orthogonalize()
{
    tree().orthogonalize();
}


template <typename T, typename Basis>
void
HTCoefficients<T, Basis>::orthogonalize_svd(std::vector<
                                            flens::DenseVector
                                            <flens::Array<T> > >& sigmas,
                                            const bool isorth)
{
    tree().orthogonalize_svd(sigmas, isorth);
}


template <typename T, typename Basis>
void
HTCoefficients<T, Basis>::truncate(const FLENS_DEFAULT_INDEXTYPE rank, bool isorth)
{
    tree().truncate(rank, isorth);
}


template <typename T, typename Basis>
void
HTCoefficients<T, Basis>::truncate(double eps)
{
    tree().truncate_hsvd(eps);
}


template <typename T, typename Basis>
T
HTCoefficients<T, Basis>::eval(const IndexD& index)
{
    assert(index.dim()==(unsigned) dim());

    typedef htucker::DimensionIndex                          IDX;

    IDX idx(dim());
    for (FLENS_DEFAULT_INDEXTYPE i=1; i<=dim(); ++i) {
        idx[i-1] = map()(index(i), i);
    }

    return httree.evaluate(idx);
}


template <typename T, typename Basis>
T
HTCoefficients<T, Basis>::eval(const IndexD& index, const FLENS_DEFAULT_INDEXTYPE vardim)
{
    assert(index.dim()==(unsigned) dim());

    htucker::DimensionIndex                 idx(dim());
    flens::DenseVector<flens::Array<FLENS_DEFAULT_INDEXTYPE> >  intindex(dim());
    for (FLENS_DEFAULT_INDEXTYPE i=1; i<=dim(); ++i) {
        intindex(i) = map()(index(i), i);
    }
    idx.setValue(intindex);

    return httree.vec_evaluate(idx, vardim);
}


template <typename T, typename Basis>
T
HTCoefficients<T, Basis>::operator()(const IndexD& index)
{
    assert(index.dim()==(unsigned) dim());
    return eval(index);
}


template <typename T, typename Basis>
T
HTCoefficients<T, Basis>::operator()(const IndexD& index,
                                     const FLENS_DEFAULT_INDEXTYPE vardim)
{
    assert(index.dim()==(unsigned) dim());
    return eval(index, vardim);
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_HTCOEFFICIENTS_TCC
