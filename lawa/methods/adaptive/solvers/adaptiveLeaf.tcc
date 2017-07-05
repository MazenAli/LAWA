#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_adaptiveLeaf_TCC
#define LAWA_METHODS_ADAPTIVE_SOLVERS_adaptiveLeaf_TCC 1

#include <cassert>
#include <iostream>
#include <cmath>

#include <lawa/methods/adaptive/righthandsides/projRhs.h>
#include <lawa/methods/adaptive/operators/operatorsd/projOperator.h>
#include <lawa/methods/adaptive/solvers/wavPcg.h>
#include <lawa/methods/adaptive/algorithms/multitreeoperations.h>

namespace lawa
{

template <typename Optype, typename Prec, typename T, typename Basis,
          typename Rhs>
unsigned
adaptiveLeaf(      Sepop<Optype>&                                A,
                   Prec&                                         P,
                   HTCoefficients<T, Basis>&                     u,
                   Coefficients<Lexicographical, T, Index1D>&    v,
                   std::vector<IndexSet<Index1D>>&               active,
                   Rhs&                                          f,
             const unsigned                                      j,
             const AdaptiveLeafParams&                           params)
{
    // Basic input check
    assert(A.dim() == (unsigned) u.dim());
    assert(A.dim() == active.size());
    assert(A.dim() == f.dim());
    assert(j>=1 && j<=A.dim());
    assert(active[j-1].size()==v.size());

    // Prepare data
    std::vector<IndexSet<Index1D> > rowvec = active;
    IndexSet<Index1D> rows                 = active[j-1];
    IndexSet<Index1D> diff                 = active[j-1];
    IndexSet<Index1D> null;
    ProjOperator<T, Optype, Basis> Ap(A, rowvec, active, u);
    Ap.setActiveDimension(j);
    Ap.computeProjection();

    // Initial residual
    Coefficients<Lexicographical, T, Index1D> b = f(null, active[j-1], j);
    Coefficients<Lexicographical, T, Index1D> r = b - Ap(v);
    T res_abs = std::sqrt(r*r);
    T res     = res_abs/b.norm(2.);

    if (params.verbose) {
        std::cout << "adaptiveLeaf: Iteration 0, r = " << res << std::endl;
    }
    if (res<=params.tol) return 0;

    // AWGM iterations
    T res_cg = 0.;
    for (unsigned k=1; k<=params.maxit; ++k) {
        // Galerkin solve
        unsigned numit = wavPcg(Ap, P, v, b, res_cg,
                                params.gamma*res_abs,
                                params.cg_maxit,
                                params.cg_verbose);

        if (params.verbose) {
            std::cout << "adaptiveLeaf: Iteration " << k
                      << ", wavPcg required " << numit
                      << " iterations to reach " << params.gamma*res
                      << std::endl;
        }

        // Extended index set
        auto oldrows = rows;
        extendMultiTree(u.basis(), diff, rows, "standard", false);
        diff = rows;
        for (const auto& it : oldrows) {
            diff.erase(it);
        }

        // Evaluate residual
        rowvec[j-1] = rows;
        b   = f(diff, rows, j);
        r   = b - Ap(v);
        P(r);
        res_abs = std::sqrt(r*r);
        res     = res_abs/b.norm(2.);

        // Check residual
        if (params.verbose) {
            std::cout << "adaptiveLeaf: Iteration " << k
                      << ", r = " << res << std::endl;
        }
        if (res<=params.tol) return k;

        // Bulk chasing
        diff = bulk(active[j-1], r, params.alpha, res_abs, res_cg, u.basis(),
                    params.bulk_verbose);

        // Update RHS and operator
        b           = f(null, active[j-1], j);
        rowvec[j-1] = active[j-1];

        // Additional information
        if (params.verbose) {
            std::cout << "adaptiveLeaf: Iteration " << k
                      << ", size of active set = "  << active[j-1].size()
                      << std::endl;
            std::cout << "adaptiveLeaf: Iteration " << k
                      << ", max active level   = "  << maxlevel(active[j-1])
                      << std::endl;
        }
    }

    std::cerr << "adaptiveLeaf: Reached max iterations " << params.maxit
              << std::endl;

    return params.maxit;
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_adaptiveLeaf_TCC
