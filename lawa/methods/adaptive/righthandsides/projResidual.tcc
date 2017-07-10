#ifndef LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_PROJRESIDUAL_TCC
#define LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_PROJRESIDUAL_TCC 1

#include <cassert>
#include <iostream>

#include <htucker/htucker.h>

namespace lawa
{

template <typename T, typename Basis, typename Optype>
ProjResidual<T, Basis, Optype>::
ProjResidual(SepCoefficients<Lexicographical, T, Index1D>& _rhsCp,
             SeparableRHSD<T, Basis>&                      _rhsInt,
             HTCoefficients<T, Basis>&                     _v,
             HTCoefficients<T, Basis>&                     _u,
             Sepop<Optype>&                                _A,
             IndexSetVec&                                  _rows,
             IndexSetVec&                                  _cols):
    rhsCp_(_rhsCp),
    rhsInt_(_rhsInt),
    v_(_v),
    u_(_u),
    A_(_A),
    rows_(_rows),
    cols_(_cols)
{
    assert(rhsCp_.dim()  == rhsInt_.dim());
    assert(rhsCp_.rank() == rhsInt_.rank());
    assert(rhsCp_.dim()  == (unsigned) v_.dim());
    assert(rhsCp_.dim()  == (unsigned) u_.dim());
    assert(rhsCp_.dim()  == A_.dim());
    assert(rhsCp_.dim()  == rows_.size());
    assert(rhsCp_.dim()  == cols_.size());
}


template <typename T, typename Basis, typename Optype>
void
ProjResidual<T, Basis, Optype>::
computeProjection(const unsigned j)
{
    assert(j>=1 && j<=dim());

    HTCoefficients<T, Basis> Au = lawa::eval(A_, u_, rows_, cols_);
    HTCoefficients<T, Basis> f(u_.dim(), u_.basis(), u_.map());
    f.tree().set_tree(u_.tree());
    set(f, rhsCp_, rows_);

    f.tree() = f.tree() - Au.tree();
    proj_    = projection(f.tree(), v_.tree(), j);
}


template <typename T, typename Basis, typename Optype>
Coefficients<Lexicographical, T, Index1D>
ProjResidual<T, Basis, Optype>::
evalLaplace(const IndexSet<Index1D>& Lambda,
            const IndexSet<Index1D>& active,
            const unsigned j)
{
    // Basic input check
    assert(A_.type()==laplace);
    assert(j>=1 && j<=dim());

    using flens::_;

    // Typedefs
    typedef SepCoefficients<Lexicographical, T, Index1D> Cp;

    // Prepare rhs
    for (unsigned i=1; i<=rhsCp_.rank(); ++i) {
        updateCoefficients(rhsCp_, i, j, rhsInt_(i, j, Lambda));
    }

    Matrix Uf  = convert(rhsCp_, active, u_, j);

    // Prepare A*u
    htucker::DimensionIndex idx(1);
    idx[0]     = j;
    Cp    _UA  = extract(u_, cols_[j-1], idx);
    Cp ret(_UA.rank(), 1);

    // Apply A*u
    for (unsigned i=1; i<=_UA.rank(); ++i) {
        TreeCoefficients1D<T> input(hashtablelength_, u_.basis().j0);
        TreeCoefficients1D<T> output(hashtablelength_, u_.basis().j0);
        typedef typename TreeCoefficients1D<T>::val_type val_type;

        // Set input
        for (const auto& lambda : _UA(i, 1)) {
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
        for (const auto& lambda : rows_[j-1]) {
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
            ret(i, 1)[lambda]  = (*it).second;
        }

        for (FLENS_DEFAULT_INDEXTYPE _j=1; _j<=JMAX; ++_j) {
            if (output.bylevel[_j].map.size()==0) break;
            for (typename CoefficientsByLevel<T>::const_it it=
                 output.bylevel[_j].map.begin();
                 it!=output.bylevel[_j].map.end(); ++it) {
                Index1D lambda(output.offset+_j, (*it).first, XWavelet);
                ret(i, 1)[lambda] = (*it).second;
            }
        }
    }

    // Convert
    int rows = maxintindhash(active, j, u_);
    Matrix UA(rows, 2*ret.rank());
    for (const auto& it : active) {
        unsigned r = u_.map()(it, j);
        for (unsigned i=1; i<=ret.rank(); ++i) {
            UA(r, i)            = _UA(i, 1)[it];
            UA(r, ret.rank()+i) = ret(i, 1)[it];
        }
    }

    // Apply projection
    Matrix f_Au(rows, UA.numCols()+Uf.numCols());
    f_Au(_, _(1, Uf.numCols())) = Uf;
    f_Au(_, _(Uf.numCols()+1, f_Au.numCols())) = UA;
    Matrix _res;
    flens::blas::mm(flens::NoTrans, flens::Trans, 1., f_Au, proj_, 0., _res);

    // Convert
    assert(_res.numCols()==1);
    return convert(_res, u_, active, j);
}


template <typename T, typename Basis, typename Optype>
Coefficients<Lexicographical, T, Index1D>
ProjResidual<T, Basis, Optype>::
eval(const IndexSet<Index1D>& Lambda,
     const IndexSet<Index1D>& active,
     const unsigned j)
{
    assert(j>=1 && j<=dim());

    Coefficients<Lexicographical, T, Index1D> ret;
    if (A_.type() != laplace) {
        std::cerr << "ProjResidual: Eval not implemented for operator type\n";
        return ret;
    }

    return evalLaplace(Lambda, active, j);
}


template <typename T, typename Basis, typename Optype>
unsigned
ProjResidual<T, Basis, Optype>::
dim() const
{
    return A_.dim();
}


template <typename T, typename Basis, typename Optype>
Coefficients<Lexicographical, T, Index1D>
ProjResidual<T, Basis, Optype>::
operator()(const IndexSet<Index1D>& Lambda,
           const IndexSet<Index1D>& active,
           const unsigned j)
{
    assert(j>=1 && j<=dim());

    return eval(Lambda, active, j);
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_PROJRESIDUAL_TCC
