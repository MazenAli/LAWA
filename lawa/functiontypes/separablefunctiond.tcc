#ifndef LAWA_FUNCTIONTYPES_SEPARABLEFUNCTIOND_TCC
#define LAWA_FUNCTIONTYPES_SEPARABLEFUNCTIOND_TCC 1

#include <cassert>

namespace lawa
{

template <typename T>
SeparableFunctionD<T>::SeparableFunctionD(const std::vector<Function<T> >& _F,
                                          const size_type          _rank,
                                          const size_type          _dim):
                     F(_F), rank_(_rank), dim_(_dim)
{
    assert(F.size()==rank_*dim_);
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
SeparableFunctionD<T>::operator()(const flens::
                                        DenseVector
                                        <flens::Array<T> >& x) const
{
    assert((unsigned) x.length()==dim_);
    typedef typename flens::DenseVector<flens::Array<T> >::IndexType IndexType;

    T res = (T) 0;
    for (size_type i=1; i<= rank_; ++i) {
        T temp = (T) 1;
        IndexType k;
        for (size_type j=1, k=x.firstIndex();
             j<=dim_; ++j, k+=x.inc()) {
            temp *= F[(j-1)*rank_+(i-1)](x(k));
        }
        res += temp;
    }

    return res;
}


template <typename T>
Function<T>&
SeparableFunctionD<T>::operator()(const size_type i, const size_type j)
{
    assert(i>=1 && i<=rank_);
    assert(j>=1 && i<=dim_);
    return F[(j-1)*rank_+(i-1)];
}


template <typename T>
const Function<T>&
SeparableFunctionD<T>::operator()(const size_type i, const size_type j) const
{
    assert(i>=1 && i<=rank_);
    assert(j>=1 && i<=dim_);
    return F[(j-1)*rank_+(i-1)];
}

} // namespace lawa

#endif // LAWA_FUNCTIONTYPES_SEPARABLEFUNCTIOND_TCC
