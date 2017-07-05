#ifndef LAWA_METHODS_ADAPTIVE_OPERATORS_OPERATORSD_PROJOPERATOR_TCC
#define LAWA_METHODS_ADAPTIVE_OPERATORS_OPERATORSD_PROJOPERATOR_TCC 1

#include <cassert>
#include <iostream>

#include <lawa/methods/adaptive/algorithms/coeffops.h>
#include <lawa/methods/adaptive/datastructures/treecoefficients1d.h>
#include <lawa/settings/enum.h>

namespace lawa
{

template <typename T, typename Optype, typename Basis>
ProjOperator<T, Optype, Basis>::
ProjOperator(Sepop<Optype>&            _A,
             IndexSetVec&              _rows,
             IndexSetVec&              _cols,
             HTCoefficients<T, Basis>& _x):
    A_(_A),
    rows_(_rows),
    cols_(_cols),
    x_(_x),
    j_(1)
{
    assert(A_.dim()==rows_.size());
    assert(A_.dim()==cols_.size());
    assert(A_.dim()==(unsigned) x_.dim());
}


template <typename T, typename Optype, typename Basis>
void
ProjOperator<T, Optype, Basis>::
setActiveDimension(const unsigned _j)
{
    assert(_j>=1 && _j<=A_.dim());
    j_ = _j;
}


template <typename T, typename Optype, typename Basis>
void
ProjOperator<T, Optype, Basis>::
computeProjection()
{
    HTCoefficients<T, Basis> Ax = lawa::eval(A_, x_, rows_, cols_);
    proj_ = projection(Ax.tree(), x_.tree(), j_);
}


template <typename T, typename Optype, typename Basis>
typename ProjOperator<T, Optype, Basis>::Matrix
ProjOperator<T, Optype, Basis>::
getProjection() const
{
    return proj_;
}


template <typename T, typename Optype, typename Basis>
Coefficients<Lexicographical, T, Index1D>
ProjOperator<T, Optype, Basis>::
eval(const Coefficients<Lexicographical, T, Index1D>& v)
{
    Coefficients<Lexicographical, T, Index1D> ret;
    if (A_.type()!=laplace) {
        std::cerr << "ProjOperator::eval: error! Operator type not implemented!\n";
        return ret;
    }

    return evalLaplace(v);
}


template <typename T, typename Optype, typename Basis>
Coefficients<Lexicographical, T, Index1D>
ProjOperator<T, Optype, Basis>::
evalLaplace(const Coefficients<Lexicographical, T, Index1D>& v)
{
    assert(A_.type()==laplace);
    assert(proj_.numRows()==1);
    assert(proj_.numCols()==2);

    Coefficients<Lexicographical, T, Index1D> ret;
    TreeCoefficients1D<T> input(hashtablelength_, x_.basis().j0);
    TreeCoefficients1D<T> output(hashtablelength_, x_.basis().j0);
    typedef typename TreeCoefficients1D<T>::val_type val_type;

    // Set input
    for (const auto& lambda : v) {
        auto _j    = lambda.first.j;
        auto _k    = lambda.first.k;
        auto xtype = lambda.first.xtype;
        if (xtype==XBSpline) {
            input.bylevel[_j-1-input.offset].
                  map.insert(val_type(_k, lambda.second));
        } else {
            input.bylevel[_j-input.offset].
                  map.insert(val_type(_k, lambda.second));
        }
    }

    // Set output
    for (const auto& lambda : rows_[j_-1]) {
        auto _j     = lambda.j;
        auto _k     = lambda.k;
        auto xtype  = lambda.xtype;
        if (xtype==XBSpline) {
            output.bylevel[_j-1-output.offset].
                  map.insert(val_type(_k, (T) 0));
        } else {
            output.bylevel[_j-output.offset].
                  map.insert(val_type(_k, (T) 0));
        }
    }

    // Apply Laplacian
    A_(1, 1).eval(input, output, "A");

    // Convert to coefficients
    for (typename CoefficientsByLevel<T>::const_it it=
            output.bylevel[0].map.begin();
            it!=output.bylevel[0].map.end(); ++it) {
        Index1D lambda(output.offset+1, (*it).first, XBSpline);
        ret[lambda]  = (*it).second;
    }

    for (FLENS_DEFAULT_INDEXTYPE i=1; i<=JMAX; ++i) {
        if (output.bylevel[i].map.size()==0) break;
        for (typename CoefficientsByLevel<T>::const_it it=
             output.bylevel[i].map.begin();
             it!=output.bylevel[i].map.end(); ++it) {
            Index1D lambda(output.offset+i, (*it).first, XWavelet);
            ret[lambda] = (*it).second;
        }
    }

    // Apply projection
    ret = proj_(1, 1)*v + proj_(1, 2)*ret;

    return ret;
}


template <typename T, typename Optype, typename Basis>
unsigned
ProjOperator<T, Optype, Basis>::
dim() const
{
    return A_.dim();
}


template <typename T, typename Optype, typename Basis>
Coefficients<Lexicographical, T, Index1D>
ProjOperator<T, Optype, Basis>::
operator()(const Coefficients<Lexicographical, T, Index1D>& v)
{
    return eval(v);
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_OPERATORS_OPERATORSD_PROJOPERATOR_TCC
