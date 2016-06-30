#include <cassert>

namespace lawa {

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
T
BasisFunction<T,Side,Domain,Cons>::operator()(T /*x*/, FLENS_DEFAULT_INDEXTYPE /*j*/, FLENS_DEFAULT_INDEXTYPE /*k*/, 
                                              unsigned short /*deriv*/) const
{
    assert(0);
    return 0.;
}

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
Support<T>
BasisFunction<T,Side,Domain,Cons>::support(FLENS_DEFAULT_INDEXTYPE /*j*/, FLENS_DEFAULT_INDEXTYPE /*k*/) const 
{
    assert(0);
    return Support<T>();
}

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
flens::DenseVector<flens::Array<T> >
BasisFunction<T,Side,Domain,Cons>::singularSupport(FLENS_DEFAULT_INDEXTYPE /*j*/, FLENS_DEFAULT_INDEXTYPE /*k*/) const 
{
    assert(0);
    return flens::DenseVector<flens::Array<T> >(); 
}

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
T
BasisFunction<T,Side,Domain,Cons>::tic(FLENS_DEFAULT_INDEXTYPE /*j*/) const
{
    assert(0);
    return 0.;
}

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
T
BasisFunction<T,Side,Domain,Cons>::getL2Norm(FLENS_DEFAULT_INDEXTYPE /*j*/, FLENS_DEFAULT_INDEXTYPE /*k*/) const
{
    assert(0);
    return 0.;
}

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
T
BasisFunction<T,Side,Domain,Cons>::getH1SemiNorm(FLENS_DEFAULT_INDEXTYPE /*j*/, FLENS_DEFAULT_INDEXTYPE /*k*/) const
{
    assert(0);
    return 0.;
}

} // namespace lawa

