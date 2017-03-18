#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_AGALS_TCC
#define LAWA_METHODS_ADAPTIVE_SOLVERS_AGALS_TCC 1

#include <algorithm>

#include <lawa/methods/adaptive/algorithms/coeffops.h>
#include <lawa/methods/adaptive/algorithms/szoneres.h>
#include <lawa/methods/adaptive/solvers/htawgm.h>

namespace lawa
{

template <typename Optype, typename Prec, typename T, typename Basis>
unsigned
agals_sym(            Engine                             *ep,
                      Sepop<Optype>&                      A,
                      Sepdiagscal<Basis>&                 S,
                      Prec&                               P,
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
    T trunc   = 1e-08;
    HTCoeff r = eval(S, F, Lambda, trunc);
    T nrmb    = nrm2(r);

    /* Initial residual */
    trunc = std::min(1e-03, params.tol*nrmb*1e-01);
    trunc   = 1e-08;
    std::cout << "trunc=" << trunc << std::endl;
    r         = eval(A, u, Lambda, Lambda);
    r.tree()  = F.tree()-r.tree();
    T unres   = nrm2(r);
    HTCoeff copy(F);
    unres    /= nrm2(copy);
    r         = eval(S, r, Lambda, trunc);
    residual  = nrm2(r)/nrmb;

    if (residual<=params.tol) {
        #ifdef VERBOSE
            std::cout << "agals_sym: Tolerance reached r = " << residual
                      << std::endl;
        #endif

        return 0;
    }

    #ifdef VERBOSE
        std::cout << "agals_sym: Iteration " << 0
                  << " r = " << residual << std::endl;
    #endif

    OptTTCoreParams p1 = params.coreopt;
    GreedyALSParams p2 = params.greedyals;
    for (unsigned k=1; k<=params.maxit; ++k) {
        unsigned greedy_it;

        /* Greedy solve */
        p1.tol    = std::min(1e-02, params.gamma*unres);
        p1.stag   = std::min(1e-04, p1.tol*1e+02);
        p2.tol    = std::min(1e-04, params.gamma*unres);
        T greedy_res;
        greedy_it = greedyALS_sym(ep, A, P, u, F, Lambda, greedy_res,
                                  params.r1update,
                                  p1,
                                  p2);
       #ifdef VERBOSE
            std::cout << "agals_sym: greedyALS required " << greedy_it
                      << " iterations to reach tolerance "
                      << greedy_res
                      << std::endl;
        #else
            (void) greedy_it;
        #endif

        /* Approximate residual */
        r         = szoneres(A, u, F, Fcp, f,
                             Lambda,
                             sweep,
                             total);
        unres     = nrm2(r);
        copy      = F;
        unres    /= nrm2(copy);
        trunc     = std::min(1e-03, residual*nrmb*1e-01);
        trunc     = 1e-08;
        std::cout << "trunc=" << trunc << std::endl;
        r         = eval(S, r, Lambda, trunc);
        residual  = nrm2(r);

        /* Bulk chasing */
        sweep = bulk(params.bulk, residual,
                     r, Lambda, sweep);

        trunc     = std::min(1e-03, residual*nrmb*1e-01);
        trunc     = 1e-08;
        r         = eval(S, F, Lambda, trunc);
        nrmb      = nrm2(r);
        residual /= nrmb;

        #ifdef VERBOSE
            std::cout << "agals_sym: Iteration " << k
                      << " r = " << residual << std::endl;
        #endif

        if (residual<=params.tol) {
            #ifdef VERBOSE
                std::cout << "agals_sym: Tolerance reached r = " << residual
                          << std::endl;
            #endif

            return k;
        }

        /* New RHS */
        restrict(F, Lambda);

        /* Extend u to new Lambda */
        rndinit(u, Lambda, 1, params.rndinit);
        //extend(u, Lambda);

        #ifdef VERBOSE
            std::cout << "agals_sym: Index set sizes\n";
            size_type size = 0;
            for (size_type j=0; j<Lambda.size(); ++j) {
                std::cout << "agals_sym: d = " << j+1
                          << " : " << Lambda[j].size() << std::endl;
                size += Lambda[j].size();
                FLENS_DEFAULT_INDEXTYPE jmax = 0;
                for (const auto& it : Lambda[j]) {
                    FLENS_DEFAULT_INDEXTYPE level = it.j;
                    if (it.xtype==XWavelet) ++level;
                    jmax = std::max(level, jmax);
                }
                std::cout << "agals_sym: jmax = " << jmax << std::endl;
            }
            std::cout << "agals_sym: Overall = " << size << std::endl;
        #endif
        u.tree().print_info();
    }

    std::cerr << "agals_sym: Max iterations reached: maxit "
              << params.maxit << " r = " << residual << std::endl;

    return params.maxit;
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_AGALS_TCC
