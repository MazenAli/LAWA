#ifndef LAWA_FUNCTIONTYPES_SEPARABLEFUNCTIOND_H
#define LAWA_FUNCTIONTYPES_SEPARABLEFUNCTIOND_H 1

#include <lawa/functiontypes/function.h>
#include <vector>
#include <cstddef>
#include <flens/flens.cxx>

namespace lawa
{

/** The RHS is of the form f = sum_{i=1^rank}\otimes_{j=1}^dim f_j^i
 *  F contains f_j^i in the order:
 *  f_1^1, f_1^2, f_1^3, ..., f_2^1, f_2^2, ..., f_dim^rank.
 */
template <typename T>
class SeparableFunctionD
{
public:
    typedef Function<T>                                     FuncType;
    typedef typename flens::DenseVector<flens::Array<T> >   DV;
    typedef typename std::vector<FuncType>::size_type       size_type;

private:
    const std::vector<FuncType>*      F;
    const size_type                   rank_;
    const size_type                   dim_;

public:
    SeparableFunctionD()                            = delete;

    SeparableFunctionD(const SeparableFunctionD&)   = delete;

    SeparableFunctionD(SeparableFunctionD&&)        = default;

    SeparableFunctionD(const std::vector<FuncType>&   _F,
                       const size_type                _rank,
                       const size_type                _dim);

    size_type
    rank() const;

    size_type
    dim() const;

    T
    eval(const DV& x) const;

    const FuncType&
    getFunction(const size_type i, const size_type j) const;

    FuncType&
    getFunction(const size_type i, const size_type j);

    T
    operator()(const DV& x) const;

    const FuncType&
    operator()(const size_type i, const size_type j) const;

    FuncType&
    operator()(const size_type i, const size_type j);
};

} // namespace lawa

#include <lawa/functiontypes/separablefunctiond.tcc>

#endif // LAWA_FUNCTIONTYPES_SEPARABLEFUNCTIOND_H
