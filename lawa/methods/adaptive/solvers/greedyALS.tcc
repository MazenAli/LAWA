#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_GREEDYALS_TCC
#define LAWA_METHODS_ADAPTIVE_SOLVERS_GREEDYALS_TCC 1

#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>

#include <htucker/htucker.h>
#include <lawa/methods/adaptive/algorithms/coeffops.h>
#include <lawa/methods/adaptive/algorithms/optTTcore.h>


namespace lawa
{

template <typename Optype, typename T, typename Basis>
unsigned
greedyALS_laplace(      Engine                             *ep,
                        Sepop<Optype>&                      A,
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
    typedef flens::GeMatrix
            <flens::FullStorage
            <T, cxxblas::ColMajor> >                                Matrix;
    typedef flens::SyMatrix
            <flens::FullStorage
            <T, cxxblas::ColMajor> >                                SyMatrix;
    typedef flens::DenseVector<
            flens::Array<FLENS_DEFAULT_INDEXTYPE> >                 IVector;

    Rank1UP_Params  p = paramsUP;
    unsigned nsweeps;
    auto bcopy = b;

    auto nrmb  = nrm2(bcopy);
    std::vector<SyMatrix>   Astiff;
    for (unsigned j=1; j<=A.dim(); ++j) {
        Astiff.push_back(assemble_projected_laplace(A, x, Lambda[j-1], j));
    }

    for (unsigned it=1; it<=params.maxit; ++it) {
        HTCoeff xold = x;
        nsweeps      = rank1update_laplace(A, Astiff, x, bcopy, Lambda, p);
        HTCoeff Ax;
        #ifdef VERBOSE
            Ax         = eval(A, x, Lambda, Lambda);
            Ax.tree()  = b.tree()-Ax.tree();
            residual   = nrm2(Ax)/nrmb;
            std::cout << "greedyALS_laplace: On update " << it
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

        Ax          = eval(A, x, Lambda, Lambda);
        Ax.tree()   = b.tree()-Ax.tree();
        residual    = nrm2(Ax)/nrmb;
        xold.tree() = x.tree()-xold.tree();
        T stag      = nrm2(xold)/nrm2(x);
        #ifdef VERBOSE
            std::cout << "greedyALS_laplace: Iteration " << it
                      << ", post core opt relative residual = "
                      << residual << std::endl;
            std::cout << "greedyALS_laplace: Iteration " << it
                      << ", post core opt        stagnation = "
                      << stag << std::endl;
        #endif





//        p.update = false;
//        nsweeps  = rank1update_laplace(A, Astiff, x, bcopy, Lambda, p);
//        #ifdef VERBOSE
//            Ax         = eval(A, x, Lambda, Lambda);
//            Ax.tree()  = b.tree()-Ax.tree();
//            residual   = nrm2(Ax)/nrmb;
//            std::cout << "greedyALS_laplace: On update " << it
//                      << " ALS required " << nsweeps
//                      << " sweeps to reach relative residual "
//                      << residual << std::endl;
//        #endif
//        p.update = true;
//
//        /* Optimize core */
//        x.tree().orthogonalize();
//        rhs = reduce_rhs(x, b);
//        B   = reduce(A, x, Lambda, Lambda);
//        htucker::extract_core(x.tree(), x0, ranks);
//        x0 = optTTcoreLaplace(ep, B, rhs, x0, ranks, paramsOpt);
//        htucker::insert_core(x.tree(), x0, ranks);
//
//        Ax          = eval(A, x, Lambda, Lambda);
//        Ax.tree()   = b.tree()-Ax.tree();
//        residual    = nrm2(Ax)/nrmb;
//        xold.tree() = x.tree()-xold.tree();
//        stag      = nrm2(xold)/nrm2(x);
//        #ifdef VERBOSE
//            std::cout << "greedyALS_laplace: Iteration " << it
//                      << ", post core opt relative residual = "
//                      << residual << std::endl;
//            std::cout << "greedyALS_laplace: Iteration " << it
//                      << ", post core opt        stagnation = "
//                      << stag << std::endl;
//        #endif

        if (residual<=params.tol || stag<=params.stag) {
            return it;
        }
    }

    std::cerr << "greedyALS_laplace: maxit " << params.maxit << " reached\n";
    return params.maxit;
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_GREEDYALS_TCC
