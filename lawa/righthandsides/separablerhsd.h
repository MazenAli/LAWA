#ifndef LAWA_RIGHTHANDSIDES_SEPARABLERHSD_H
#define LAWA_RIGHTHANDSIDES_SEPARABLERHSD_H 1

#include <lawa/functiontypes/separablefunctiond.h>
#include <vector>
#include <flens/flens.cxx>
#include <lawa/integrals/integrals.h>
#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/datastructures/index.h>

namespace lawa
{

template <typename T, typename Basis>
class SeparableRHSD
{
public:
    typedef IntegralF<Gauss, Basis>                             IntegralType;
    typedef SeparableFunctionD<T>                               FuncType;
    typedef typename flens::GeMatrix<flens::FullStorage
                                    <T, cxxblas::ColMajor> >    GeMat;
    typedef typename flens::GeMatrix<flens::FullStorage
                                    <int, cxxblas::ColMajor> >  IntMat;
    typedef typename std::vector<GeMat>                         Matvec;
    typedef typename std::vector<IntegralType>                  Intvec;

    typedef typename Matvec::size_type                          size_type_Mat;
    typedef typename Intvec::size_type                          size_type;
    typedef Coefficients<Lexicographical, T, Index1D>           Coeff1D;

private:
    const size_type                 rank_;
    const size_type                 dim_;
    const Basis*                    basis;
    const FuncType*                 F;
    const Matvec*                   deltas;
    const IntMat*                   derivs;
    Intvec                          integral;

public:
    SeparableRHSD()                     = delete;

    SeparableRHSD(const SeparableRHSD&) = delete;

    SeparableRHSD(SeparableRHSD&&)      = default;

    SeparableRHSD(const Basis& _basis, const FuncType& _F,
                  const int order = 4);

    SeparableRHSD(const Basis& _basis,   const FuncType& _F,
                  const Matvec& _deltas, const IntMat& _derivs,
                  const int order = 4);

    size_type
    rank() const;

    size_type
    dim() const;

    const IntegralType&
    getInt(const size_type i, const size_type j) const;

    IntegralType&
    getInt(const size_type i, const size_type j);

    const GeMat&
    getDeltas(const size_type i, const size_type j) const;

    T
    eval(const size_type i, const size_type j, const Index1D& index) const;

    Coeff1D
    eval(const size_type i, const size_type j, const IndexSet<Index1D>&
         indexset) const;

    T
    eval(const IndexD& index) const;

    T
    operator()(const size_type i, const size_type j,
               const Index1D& index) const;

    Coeff1D
    operator()(const size_type i, const size_type j,
               const IndexSet<Index1D>& index) const;

    T
    operator()(const IndexD& index) const;
};

} // namespace lawa

#include <lawa/righthandsides/separablerhsd.tcc>

#endif // LAWA_RIGHTHANDSIDES_SEPARABLERHSD_H
