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
    typedef std::vector<Coeff1D>            mapType;

    const size_type                 rank_;
    const size_type                 dim_;
    const Basis*                    basis;
    const FuncType*                 F;
    const Matvec*                   deltas;
    const IntMat*                   derivs;
    Intvec                          integral;
    mapType                         data;

public:
    SeparableRHSD()                     = delete;

    SeparableRHSD(const SeparableRHSD&) = delete;

    SeparableRHSD(SeparableRHSD&&)      = default;

    SeparableRHSD(const Basis& _basis, const FuncType& _F,
                  const FLENS_DEFAULT_INDEXTYPE order = 4);

    SeparableRHSD(const Basis& _basis,   const FuncType& _F,
                  const Matvec& _deltas, const IntMat& _derivs,
                  const FLENS_DEFAULT_INDEXTYPE order = 4);

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
    eval(const size_type i, const size_type j, const Index1D& index);

    Coeff1D
    eval(const size_type i, const size_type j, const IndexSet<Index1D>&
         indexset);

    T
    eval(const IndexD& index);

    T
    operator()(const size_type i, const size_type j,
               const Index1D& index);

    Coeff1D
    operator()(const size_type i, const size_type j,
               const IndexSet<Index1D>& index);

    T
    operator()(const IndexD& index);
};

} // namespace lawa

#include <lawa/righthandsides/separablerhsd.tcc>

#endif // LAWA_RIGHTHANDSIDES_SEPARABLERHSD_H
