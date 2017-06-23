#ifndef LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_PROJRHS_TCC
#define LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_PROJRHS_TCC 1

#include <cassert>

namespace lawa
{

template <typename T, typename Basis>
ProjRhs<T, Basis>::
ProjRhs(      SepCoefficients<Lexicographical, T, Index1D>& _rhsCp,
        const SeparableRHSD<T, Basis>&                      _rhsInt,
              HTCoefficients<T, Basis>&                     _tree):
    rhsCp_(_rhsCp),
    rhsInt_(_rhsInt),
    tree_(_tree)
{}


template <typename T, typename Basis>
void
ProjRhs<T, Basis>::
computeProjection(const unsigned j)
{
    assert(j>=1 && j<=rhsCp_.dim());

    HTCoefficients<T, Basis> b(tree_.dim(), tree_.basis(), tree_.map());
    set(b, rhsCp_);
    proj_ = projection(b.tree(), tree_.tree(), j);
}


template <typename T, typename Basis>
Coefficients<Lexicographical, T, Index1D>
ProjRhs<T, Basis>::
eval(const IndexSet<Index1D>& Lambda,
     const IndexSet<Index1D>& active,
     const unsigned           j)
{
    assert(j>=1 && j<=rhsCp_.dim());

    for (unsigned i=1; i<=rhsCp_.rank(); ++i) {
        addCoefficients(rhsCp_, i, j, rhsInt_(i, j, Lambda));
    }

    Matrix U = convert(rhsCp_, tree_, j);

    Matrix ret;
    flens::blas::mm(flens::NoTrans, flens::Trans, 1., U, proj_, 0., ret);

    Coefficients<Lexicographical, T, Index1D>
    out = convert(ret, tree_, active, j);

    return out;
}


template <typename T, typename Basis>
typename ProjRhs<T, Basis>::Matrix
ProjRhs<T, Basis>::
getProj() const
{
    return proj_;
}


template <typename T, typename Basis>
Coefficients<Lexicographical, T, Index1D>
ProjRhs<T, Basis>::
operator()(const IndexSet<Index1D>& Lambda,
           const IndexSet<Index1D>& active,
           const unsigned           j)
{
    assert(j>=1 && j<=rhsCp_.dim());

    return eval(Lambda, active, j);
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_PROJRHS_TCC
