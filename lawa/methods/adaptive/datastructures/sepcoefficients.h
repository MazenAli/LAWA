#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_SEPCOEFFICIENTS_H
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_SEPCOEFFICIENTS_H 1

#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/righthandsides/separablerhsd.h>
#include <vector>
#include <iostream>

namespace lawa
{

template <SortingCriterion S, typename T, typename Index>
class SepCoefficients
{
public:
    typedef Coefficients<S, T, Index>       Coeff;
    typedef typename std::vector<Coeff>     CoeffVec;
    typedef typename CoeffVec::size_type    size_type;

private:
    size_type     rank_;
    size_type     dim_;
    CoeffVec      coeffs;

public:
    SepCoefficients();

    SepCoefficients(const SepCoefficients&) = default;

    SepCoefficients(SepCoefficients&&)      = default;

    SepCoefficients(const size_type _rank, const size_type _dim);

    SepCoefficients(const CoeffVec& _coeffs,
                    const size_type _rank, const size_type _dim);

    void
    resize(const size_type _rank, const size_type _dim);

    size_type
    rank() const;

    size_type
    dim() const;

    const CoeffVec&
    getCoefficients() const;

    CoeffVec&
    getCoefficients();

    const Coeff&
    getCoefficients(const size_type i, const size_type j) const;

    Coeff&
    getCoefficients(const size_type i, const size_type j);

    T
    eval(const size_type i, const size_type j,
         const Index1D& index) const;

    T
    eval(const IndexD& index) const;

    const Coeff&
    operator()(const size_type i, const size_type j) const;

    Coeff&
    operator()(const size_type i, const size_type j);

    T
    operator()(const size_type i, const size_type j,
               const Index1D& index) const;

    T
    operator()(const IndexD& index) const;

    SepCoefficients<S, T, Index>&
    operator=(const SepCoefficients<S, T, Index>& copy);
};

} // namespace lawa

#include <lawa/methods/adaptive/datastructures/sepcoefficients.tcc>

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_SEPCOEFFICIENTS_H
