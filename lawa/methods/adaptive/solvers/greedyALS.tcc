#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_GREEDYALS_TCC
#define LAWA_METHODS_ADAPTIVE_SOLVERS_GREEDYALS_TCC 1

#include <iostream>
#include <lawa/methods/adaptive/algorithms/coeffops.h>
#include <lawa/methods/adaptive/algorithms/optTTcore.h>
#include <cassert>
#include <htucker/htucker.h>

namespace lawa
{

template <typename Optype, typename Prec, typename T, typename Basis>
unsigned
greedyALS_sym(        Engine                             *ep,
                      Sepop<Optype>&                      A,
                      Prec&                               P,
                      HTCoefficients<T, Basis>&           x,
                const HTCoefficients<T, Basis>&           b,
                const std::vector<IndexSet<Index1D> >&    Lambda,
                      T&                                  residual,
                const Rank1UP_Params&                     paramsUP,
                const OptTTCoreParams&                    paramsOpt,
                const GreedyALSParams&                    params)
{
    assert(A.dim()==(unsigned) x.dim());
    assert(A.dim()==(unsigned) b.dim());
    assert(A.dim()==Lambda.size());

    typedef HTCoefficients<T, Basis>                                HTCoeff;
    typedef flens::GeMatrix<flens::FullStorage
            <T, cxxblas::ColMajor> >                                Matrix;
    typedef flens::DenseVector<
            flens::Array<FLENS_DEFAULT_INDEXTYPE> >                 IVector;

    unsigned it;
    Rank1UP_Params  p = paramsUP;
    unsigned nsweeps;
    auto bcopy = b;
    auto nrmb  = nrm2(bcopy);
    for (it=1; it<=params.maxIt; ++it) {
        nsweeps = rank1update_sym(A, P, x, bcopy, Lambda, p);
        HTCoeff Ax;
        #ifdef VERBOSE
            Ax         = eval(A, x, Lambda, Lambda);
            Ax.tree()  = b.tree()-Ax.tree();
            residual   = nrm2(Ax)/nrmb;
            std::cout << "greedyALS_sym: On update " << it
                      << " ALS required " << nsweeps
                      << " sweeps to reach relative residual "
                      << residual << std::endl;
        #endif
        p.update = true;

        /* Optimize core */
        x.tree().orthogonalize();
        auto rhs = reduce_rhs(x, b);
        auto B   = reduce(A, x, Lambda, Lambda);
        std::vector<Matrix> x0(rhs.size()-1);
        IVector ranks(x0.size());
        htucker::extract_core(x.tree(), x0, ranks);
        x0 = optTTcoreLaplace(ep, B, rhs, x0, ranks, paramsOpt);
        htucker::insert_core(x.tree(), x0, ranks);

        Ax         = eval(A, x, Lambda, Lambda);
        Ax.tree()  = b.tree()-Ax.tree();
        residual   = nrm2(Ax)/nrmb;
        #ifdef VERBOSE
            std::cout << "greedyALS_sym: Post core opt relative residual = "
                      << residual << std::endl;
        #endif

        if (residual<=params.tol) {
            return it;
        }
    }

    std::cerr << "greedyALS_sym: maxit " << params.maxIt << " reached\n";
    return it;
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_GREEDYALS_TCC
