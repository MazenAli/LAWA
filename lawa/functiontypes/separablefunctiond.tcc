#ifndef LAWA_FUNCTIONTYPES_SEPARABLEFUNCTIOND_TCC
#define LAWA_FUNCTIONTYPES_SEPARABLEFUNCTIOND_TCC 1

#include <cassert>

namespace lawa
{

template <typename T>
SeparableFunctionD<T>::SeparableFunctionD(const std::vector<FuncType>& _F,
                                          const size_type          _rank,
                                          const size_type          _dim):
    F(&_F),
    rank_(_rank),
    dim_(_dim)
{
    assert(F->size()==rank()*dim());
}


template <typename T>
typename SeparableFunctionD<T>::size_type
SeparableFunctionD<T>::rank() const
{
    return rank_;
}


template <typename T>
typename SeparableFunctionD<T>::size_type
SeparableFunctionD<T>::dim() const
{
    return dim_;
}


template <typename T>
T
SeparableFunctionD<T>::eval(const DV& x) const
{
    assert((size_type) x.length()==dim());
    typedef typename flens::DenseVector<flens::Array<T> >::IndexType IndexType;

    T res = (T) 0;
    for (size_type i=1; i<= rank(); ++i) {
        T temp = (T) 1;
        IndexType k;
        for (size_type j=1, k=x.firstIndex();
             j<=dim(); ++j, k+=x.inc()) {
            temp *= getFunction(i, j)(x(k));
        }
        res += temp;
    }

    return res;
}


template <typename T>
const typename SeparableFunctionD<T>::FuncType&
SeparableFunctionD<T>::getFunction(const size_type i, const size_type j) const
{
    assert(i>=1 && i<=rank());
    assert(j>=1 && i<=dim());
    return (*F)[(j-1)*rank()+(i-1)];
}


template <typename T>
typename SeparableFunctionD<T>::FuncType&
SeparableFunctionD<T>::getFunction(const size_type i, const size_type j)
{
    assert(i>=1 && i<=rank());
    assert(j>=1 && i<=dim());
    return (*F)[(j-1)*rank()+(i-1)];
}


template <typename T>
T
SeparableFunctionD<T>::operator()(const DV& x) const
{
    return eval(x);
}


template <typename T>
const typename SeparableFunctionD<T>::FuncType&
SeparableFunctionD<T>::operator()(const size_type i, const size_type j) const
{
    assert(i>=1 && i<=rank());
    assert(j>=1 && i<=dim());
    return getFunction(i, j);
}


template <typename T>
typename SeparableFunctionD<T>::FuncType&
SeparableFunctionD<T>::operator()(const size_type i, const size_type j)
{
    assert(i>=1 && i<=rank());
    assert(j>=1 && i<=dim());
    return getFunction(i, j);
}

} // namespace lawa

#endif // LAWA_FUNCTIONTYPES_SEPARABLEFUNCTIOND_TCC
