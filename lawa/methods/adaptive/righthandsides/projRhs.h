#ifndef LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_PROJRHS_H
#define LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_PROJRHS_H 1

#include <flens/flens.cxx>

#include <lawa/methods/adaptive/datastructures/sepcoefficients.h>
#include <lawa/righthandsides/separablerhsd.h>
#include <lawa/methods/adaptive/datastructures/htcoefficients.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/index.h>

namespace lawa
{

// A wrapper class for the RHS for adaptive ALS routines
template <typename T, typename Basis>
class ProjRhs
{

public:
    typedef typename flens::GeMatrix<
                     flens::FullStorage<T, flens::ColMajor> > Matrix;
    typedef std::vector<IndexSet<Index1D> >                   IndexSetVec;

    ProjRhs(      SepCoefficients<Lexicographical, T, Index1D>& _rhsCp,
                  IndexSetVec&                                  _active,
                  SeparableRHSD<T, Basis>&                      _rhsInt,
                  HTCoefficients<T, Basis>&                     _tree);

    void
    computeProjection(const unsigned j);

    Coefficients<Lexicographical, T, Index1D>
    eval(const IndexSet<Index1D>& Lambda,
         const IndexSet<Index1D>& active,
         const unsigned           j);

    Matrix
    getProj() const;

    unsigned
    dim() const;

    Coefficients<Lexicographical, T, Index1D>
    operator()(const IndexSet<Index1D>& Lambda,
               const IndexSet<Index1D>& active,
               const unsigned           j);

private:
    SepCoefficients<Lexicographical, T, Index1D>& rhsCp_;
    IndexSetVec&                                  active_;
    SeparableRHSD<T, Basis>&                      rhsInt_;
    HTCoefficients<T, Basis>&                     tree_;
    Matrix                                        proj_;

}; // class ProjRhs

} // namespace lawa

#include <lawa/methods/adaptive/righthandsides/projRhs.tcc>

#endif // LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_PROJRHS_H
