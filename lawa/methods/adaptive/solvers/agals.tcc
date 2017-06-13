#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_AGALS_TCC
#define LAWA_METHODS_ADAPTIVE_SOLVERS_AGALS_TCC 1

#include <algorithm>

#include <lawa/methods/adaptive/algorithms/coeffops.h>
#include <lawa/methods/adaptive/algorithms/szoneres.h>
#include <lawa/methods/adaptive/solvers/htawgm.h>

#include <extensions/general/vmemusage.h>

namespace lawa
{

template <typename Optype, typename T, typename Basis>
unsigned
agals_laplace(      Engine                             *ep,
                    Sepop<Optype>&                      A,
                    Sepdiagscal<Basis>&                 S,
                    HTCoefficients<T, Basis>&           u,
              const SeparableRHSD<T, Basis>&            f,
                    std::vector<IndexSet<Index1D> >&    Lambda,
                    T&                                  residual,
              const AgALSParams&                        params)
{
    assert(A.dim()==S.dim());
    assert(A.dim()==(unsigned) u.dim());
    assert(A.dim()==Lambda.size());

    typedef HTCoefficients<T, Basis>                    HTCoeff;
    typedef std::vector<IndexSet<Index1D> >::size_type  size_type;

    std::vector<IndexSet<Index1D> >              sweep(Lambda);
    std::vector<IndexSet<Index1D> >              total(Lambda);

    SepCoefficients<Lexicographical, T, Index1D> Fcp(f.rank(), f.dim());
    HTCoeff                                      F(u.dim(), u.basis(), u.map());
    F.tree().set_tree(u.tree());

    genCoefficients(Fcp, f, Lambda);
    set(F, Fcp);
    T trunc   = 1e-04;
    HTCoeff r = applyScaleTT(S, F, Lambda, trunc);
    T nrmb    = nrm2(r);

    /* Initial residual */
    r         = eval(A, u, Lambda, Lambda);
    r.tree()  = F.tree()-r.tree();
    T unres   = nrm2(r);
    HTCoeff copy(F);
    unres    /= nrm2(copy);
    r         = applyScaleTT(S, r, Lambda, trunc);
    residual  = nrm2(r)/nrmb;

    if (residual<=params.tol) {
        #ifdef VERBOSE
            std::cout << "agals_laplace: Tolerance reached r = " << residual
                      << std::endl;
        #endif

        return 0;
    }

    #ifdef VERBOSE
        std::cout << "agals_laplace: Iteration " << 0
                  << " r = " << residual << std::endl;
    #endif

    Rank1UP_Params  p0 = params.r1update;
    OptTTCoreParams p1 = params.coreopt;
    GreedyALSParams p2 = params.greedyals;
    for (unsigned k=1; k<=params.maxit; ++k) {
        unsigned greedy_it;
        size_type size = 0;

        /* Greedy solve */
        p1.tol    = std::min(1e-06, params.gamma*residual);
        p1.stag   = 1e-03;
        T greedy_res;
        greedy_it = greedyALS_laplace(ep, A, u, F, Lambda, greedy_res,
                                      p0,
                                      p1,
                                      p2);

       #ifdef VERBOSE
            std::cout << "agals_laplace: greedyALS required " << greedy_it
                      << " iterations to reach tolerance "
                      << greedy_res
                      << std::endl;
        #else
            (void) greedy_it;
        #endif

        /* Approximate residual */
        r = szoneres(A, u, F, Fcp, f,
                     Lambda,
                     sweep,
                     total);
        trunc     = std::min(1e-03, residual*nrmb*1e-01);
        if (trunc<1e-10) trunc = 1e-10;
        S.set_nu(residual);
        r         = applyScaleTT(S, r, total, trunc);
        residual  = nrm2(r);

        /* Bulk chasing */
        sweep = bulk(params.bulk, residual,
                     r, Lambda, sweep);

        trunc     = std::min(1e-03, residual*nrmb*1e-01);
        restrict(F, Lambda);
        r         = applyScaleTT(S, F, Lambda, trunc);
        nrmb      = nrm2(r);
        residual /= nrmb;

        #ifdef VERBOSE
            std::cout << "agals_laplace: Iteration " << k
                      << " r = " << residual << std::endl;
        #endif

        if (residual<=params.tol) {
            #ifdef VERBOSE
                std::cout << "agals_laplace: Tolerance reached r = " << residual
                          << std::endl;
            #endif

            return k;
        }

        /* Extend u to new Lambda */
        extend(u, Lambda);

        /* Post smoothing */
//        if (k>1) {
//            std::cout << "agals_laplace: Post extend smoothing...\n";
//
//            p0.update  = false;
//            (void) greedyALS_laplace(ep, A, u, F, Lambda, greedy_res,
//                                     p0,
//                                     p1,
//                                     p2);
//
//        }
        p0.update = true;

        #ifdef VERBOSE
            std::cout << "agals_laplace: Current max rank of solution "
                      << u.tree().max_rank() << std::endl;
            std::cout << "agals_laplace: Index set sizes\n";
            size = 0;
            for (size_type j=0; j<Lambda.size(); ++j) {
                std::cout << "agals_laplace: d = " << j+1
                          << " : " << Lambda[j].size() << std::endl;
                size += Lambda[j].size();
                FLENS_DEFAULT_INDEXTYPE jmax = 0;
                for (const auto& it : Lambda[j]) {
                    FLENS_DEFAULT_INDEXTYPE level = it.j;
                    if (it.xtype==XWavelet) ++level;
                    jmax = std::max(level, jmax);
                }
                std::cout << "agals_laplace: jmax = " << jmax << std::endl;
            }
            std::cout << "agals_laplace: Overall = " << size << std::endl;
        #endif
    }

    std::cerr << "agals_laplace: Max iterations reached: maxit "
              << params.maxit << " r = " << residual << std::endl;

    return params.maxit;
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_AGALS_TCC
