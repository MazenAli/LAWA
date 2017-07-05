#ifndef LAWA_METHODS_ADAPTIVE_OPERATORS_OPERATORSD_PROJOPERATOR_H
#define LAWA_METHODS_ADAPTIVE_OPERATORS_OPERATORSD_PROJOPERATOR_H 1

#include <cstddef>
#include <vector>

#include <flens/flens.cxx>

#include <lawa/methods/adaptive/operators/operatorsd/sepop.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/datastructures/htcoefficients.h>

namespace lawa
{

// A wrapper class for the operator for ALS routines
template <typename T, typename Optype, typename Basis>
class ProjOperator
{
public:
    typedef typename flens::GeMatrix<
                     flens::FullStorage<T, flens::ColMajor> >  Matrix;
    typedef std::vector<IndexSet<Index1D> >                    IndexSetVec;

    ProjOperator(Sepop<Optype>&            _A,
                 IndexSetVec&              _rows,
                 IndexSetVec&              _cols,
                 HTCoefficients<T, Basis>& _x);

    void
    setActiveDimension(const unsigned _j);

    void
    computeProjection();

    Matrix
    getProjection() const;

    Coefficients<Lexicographical, T, Index1D>
    eval(const Coefficients<Lexicographical, T, Index1D>& v);

    Coefficients<Lexicographical, T, Index1D>
    evalLaplace(const Coefficients<Lexicographical, T, Index1D>& v);

    unsigned
    dim() const;

    Coefficients<Lexicographical, T, Index1D>
    operator()(const Coefficients<Lexicographical, T, Index1D>& v);

private:
    Sepop<Optype>&            A_;
    IndexSetVec&              rows_;
    IndexSetVec&              cols_;
    HTCoefficients<T, Basis>& x_;
    unsigned                  j_;
    Matrix                    proj_;
    std::size_t               hashtablelength_ = 193;
};

} // namespace lawa

#include <lawa/methods/adaptive/operators/operatorsd/projOperator.tcc>

#endif // LAWA_METHODS_ADAPTIVE_OPERATORS_OPERATORSD_PROJOPERATOR_H
