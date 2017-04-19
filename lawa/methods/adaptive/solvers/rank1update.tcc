#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_RANK1UPDATE_TCC
#define LAWA_METHODS_ADAPTIVE_SOLVERS_RANK1UPDATE_TCC 1

#include <iostream>
#include <lawa/methods/adaptive/algorithms/coeffops.h>
#include <cassert>

#include <htucker/htucker.h>

namespace lawa
{

template <typename Optype, typename Prec, typename T, typename Basis>
unsigned
rank1update_sym(      Sepop<Optype>&                      A,
                      Prec&                               P,
                      HTCoefficients<T, Basis>&           x,
                const HTCoefficients<T, Basis>&           b,
                const std::vector<IndexSet<Index1D> >&    Lambda,
                const flens::GeMatrix
                      <flens::FullStorage
                      <T, cxxblas::ColMajor> >&           scales,
                const Rank1UP_Params&                     params)
{
    assert(A.dim()==(unsigned) x.dim());
    assert(A.dim()==(unsigned) b.dim());
    assert(A.dim()==Lambda.size());

    typedef HTCoefficients<T, Basis>                                HTCoeff;

    HTCoeff bup(x.dim(), x.basis(), x.map());
    bup.tree().set_tree(x.tree());
    if (params.update) {
        HTCoeff Ax = eval(A, x, Lambda, Lambda);
        bup.tree() = b.tree()-Ax.tree();
    }


    unsigned it;
    T residual;
    if (params.update) {
        HTCoeff xup(x.dim(), x.basis(), x.map());
        xup.tree().set_tree(x.tree());
        rndinit(xup, Lambda, 1, 1e-03);
        it = precrank1als_sym(A, P, xup, bup, Lambda, scales, residual,
                              params.check_res,
                              params.orthog,
                              params.sw,
                              params.balance,
                              params.tol_als,
                              params.max_sweep,
                              params.tol_cg,
                              params.maxit_cg,
                              params.verbose);
        x.tree() = x.tree()+xup.tree();
    } else {
        it = precrank1als_sym(A, P, x, b, Lambda, scales, residual,
                              params.check_res,
                              params.orthog,
                              params.sw,
                              params.balance,
                              params.tol_als,
                              params.max_sweep,
                              params.tol_cg,
                              params.maxit_cg,
                              params.verbose);
    }

    return it;
}


template <typename Optype, typename T, typename Basis>
unsigned
rank1update_laplace(      Sepop<Optype>&                      A,
                    const std::vector<flens::SyMatrix
                          <flens::FullStorage
                          <T, cxxblas::ColMajor> > >&         Astiff,
                          HTCoefficients<T, Basis>&           x,
                    const HTCoefficients<T, Basis>&           b,
                    const std::vector<IndexSet<Index1D> >&    Lambda,
                    const Rank1UP_Params&                     params)
{
    assert(A.dim()==(unsigned) x.dim());
    assert(A.dim()==(unsigned) b.dim());
    assert(A.dim()==Lambda.size());

    typedef HTCoefficients<T, Basis>                                HTCoeff;

    HTCoeff bup(x.dim(), x.basis(), x.map());
    bup.tree().set_tree(x.tree());
    if (params.update) {
        HTCoeff Ax = eval(A, x, Lambda, Lambda);
        bup.tree() = b.tree()-Ax.tree();
    }


    unsigned it;
    if (params.update) {
        HTCoeff xup(x.dim(), x.basis(), x.map());
        xup.tree().set_tree(x.tree());
        rndinit(xup, Lambda, 1, 1e-03);
        it = rank1als_laplace(A, Astiff, xup, bup, Lambda,
                              params.check_res,
                              params.orthog,
                              params.tol_als,
                              params.max_sweep,
                              params.verbose);
        x.tree() = x.tree()+xup.tree();
    } else {
        it = rank1als_laplace(A, Astiff, x, b, Lambda,
                              params.check_res,
                              params.orthog,
                              params.tol_als,
                              params.max_sweep,
                              params.verbose);
    }

    return it;
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_RANK1UPDATE_TCC
