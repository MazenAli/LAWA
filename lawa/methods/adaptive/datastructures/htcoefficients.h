#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_HTCOEFFICIENTS_H
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_HTCOEFFICIENTS_H 1

#include <htucker/htucker.h>
#include <cstddef>
#include <lawa/methods/adaptive/datastructures/sepcoefficients.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/mapwavind.h>
#include <flens/flens.cxx>
#include <vector>

namespace lawa
{

template <typename T, typename Basis>
class HTCoefficients
{
public:
    typedef typename htucker::HTuckerTree<T>    HTTree;
    typedef Mapwavind<Index1D>                  Maptype;

private:
    HTTree              httree;
    const Basis*        basis_;
    Maptype*            map_;

public:
    HTCoefficients();

    HTCoefficients(const HTCoefficients&)   = default;

    HTCoefficients(HTCoefficients&&)        = default;

    HTCoefficients(const FLENS_DEFAULT_INDEXTYPE d, const Basis& _basis,
                   Maptype& _map);

    HTCoefficients(const FLENS_DEFAULT_INDEXTYPE d, const double split,
                   const Basis& _basis, Maptype& _map);

    FLENS_DEFAULT_INDEXTYPE
    dim() const;

    void
    print() const;

    const HTTree&
    tree() const;

    HTTree&
    tree();

    const Basis&
    basis() const;

    const Maptype&
    map() const;

    Maptype&
    map();

    void
    orthogonalize();

    void
    orthogonalize_svd(std::vector<
                      flens::DenseVector<flens::Array<T> > >& sigmas,
                      const bool isorth = false);

    void
    truncate(const FLENS_DEFAULT_INDEXTYPE rank, bool isorth = false);

    void
    truncate(double eps);

    T
    eval(const IndexD& index);

    T
    eval(const IndexD& index, const FLENS_DEFAULT_INDEXTYPE vardim);

    T
    operator()(const IndexD& index);

    T
    operator()(const IndexD& index, const FLENS_DEFAULT_INDEXTYPE vardim);

    HTCoefficients<T, Basis>&
    operator=(const HTCoefficients<T, Basis>& copy);
};

} // namespace lawa

#include <lawa/methods/adaptive/datastructures/htcoefficients.tcc>

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_HTCOEFFICIENTS_H
