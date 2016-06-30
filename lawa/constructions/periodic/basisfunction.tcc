#include <cassert>

namespace lawa {

template <typename T, FunctionSide Side, Construction Cons>
T
BasisFunction<T,Side,Periodic,Cons>::operator()(T /*x*/, FLENS_DEFAULT_INDEXTYPE /*j*/, FLENS_DEFAULT_INDEXTYPE /*k*/, 
                                                unsigned short /*deriv*/) const
{
    assert(0);
    return 0.;
}

template <typename T, FunctionSide Side, Construction Cons>
PeriodicSupport<T>
BasisFunction<T,Side,Periodic,Cons>::support(FLENS_DEFAULT_INDEXTYPE /*j*/, FLENS_DEFAULT_INDEXTYPE /*k*/) const 
{
    assert(0);
    return PeriodicSupport<T>();
}

template <typename T, FunctionSide Side, Construction Cons>
flens::DenseVector<flens::Array<T> >
BasisFunction<T,Side,Periodic,Cons>::singularSupport(FLENS_DEFAULT_INDEXTYPE /*j*/, FLENS_DEFAULT_INDEXTYPE /*k*/) const 
{
    assert(0);
    return flens::DenseVector<flens::Array<T> >(); 
}

template <typename T, FunctionSide Side, Construction Cons>
T
BasisFunction<T,Side,Periodic,Cons>::tic(FLENS_DEFAULT_INDEXTYPE /*j*/) const
{
    assert(0);
    return 0.;
}

} // namespace lawa

