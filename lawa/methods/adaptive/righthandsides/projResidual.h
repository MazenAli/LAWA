#ifndef LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_PROJRESIDUAL_H
#define LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_PROJRESIDUAL_H 1

#include <vector>

#include <flens/flens.cxx>

#include <lawa/methods/adaptive/datastructures/sepcoefficients.h>
#include <lawa/righthandsides/separablerhsd.h>
#include <lawa/methods/adaptive/datastructures/htcoefficients.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/index.h>

namespace lawa
{

// A wrapper class for the RHS for rank 1 ALS updates
template <typename T, typename Basis, typename Optype>
class ProjResidual
{

public:
    typedef typename flens::GeMatrix<
                     flens::FullStorage<T, flens::ColMajor> > Matrix;
    typedef std::vector<IndexSet<Index1D> >                   IndexSetVec;

    ProjResidual(SepCoefficients<Lexicographical, T, Index1D>& _rhsCp,
                 SeparableRHSD<T, Basis>&                      _rhsInt,
                 HTCoefficients<T, Basis>&                     _v,
                 HTCoefficients<T, Basis>&                     _u,
                 Sepop<Optype>&                                _A,
                 IndexSetVec&                                  _rows,
                 IndexSetVec&                                  _cols);

    void
    computeProjection(const unsigned j);

    Coefficients<Lexicographical, T, Index1D>
    evalLaplace(const IndexSet<Index1D>& Lambda,
                const IndexSet<Index1D>& active,
                const unsigned j);


    Coefficients<Lexicographical, T, Index1D>
    eval(const IndexSet<Index1D>& Lambda,
         const IndexSet<Index1D>& active,
         const unsigned j);

    unsigned
    dim() const;

    Coefficients<Lexicographical, T, Index1D>
    operator()(const IndexSet<Index1D>& Lambda,
               const IndexSet<Index1D>& active,
               const unsigned j);

private:
    SepCoefficients<Lexicographical, T, Index1D>& rhsCp_;
    SeparableRHSD<T, Basis>&                      rhsInt_;
    HTCoefficients<T, Basis>&                     v_;
    HTCoefficients<T, Basis>&                     u_;
    Sepop<Optype>&                                A_;
    IndexSetVec&                                  rows_;
    IndexSetVec&                                  cols_;
    std::size_t                                   hashtablelength_ = 193;
    Matrix                                        proj_;

}; // clas ProjResidual

} // namespace lawa

#include <lawa/methods/adaptive/righthandsides/projResidual.tcc>

#endif // LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_PROJRESIDUAL_H
