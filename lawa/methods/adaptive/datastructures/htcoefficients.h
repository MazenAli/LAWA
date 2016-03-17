#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_HTCOEFFICIENTS_H
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_HTCOEFFICIENTS_H 1

#include <htucker/htucker.h>
#include <lawa/methods/adaptive/datastructures/sepcoefficients.h>
#include <lawa/methods/adaptive/datastructures/index.h>

namespace lawa
{

template <typename T, typename Basis>
class HTCoefficients
{
public:
    typedef typename htucker::HTuckerTree<T>    HTTree;

private:
    HTTree          httree;
    const Basis*    basis_;

public:
    HTCoefficients();

    HTCoefficients(const HTCoefficients&)   = default;

    HTCoefficients(HTCoefficients&&)        = default;

    HTCoefficients(const int d, const Basis& _basis);

    HTCoefficients(const int d, const double split, const Basis& _basis);

    int
    dim() const;

    void
    print() const;

    const HTTree&
    tree() const;

    HTTree&
    tree();

    const Basis&
    basis() const;

    void
    orthogonolize();

    void
    truncate(const int rank, bool isorth = false);

    void
    truncate(double eps, bool isorth = false);

    T
    eval(const IndexD& index) const;

    T
    eval(const IndexD& index, const int vardim) const;

    T
    operator()(const IndexD& index) const;

    T
    operator()(const IndexD& index, const int vardim) const;

    HTCoefficients<T, Basis>&
    operator=(const HTCoefficients<T, Basis>& copy);
};

} // namespace lawa

#include <lawa/methods/adaptive/datastructures/htcoefficients.tcc>

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_HTCOEFFICIENTS_H
