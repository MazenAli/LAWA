#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_RANK1ADAPTIVEALS_TCC
#define LAWA_METHODS_ADAPTIVE_SOLVERS_RANK1ADAPTIVEALS_TCC 1

#include <cassert>
#include <iostream>

#include <lawa/methods/adaptive/solvers/adaptiveLeaf.h>
#include <lawa/methods/adaptive/algorithms/coeffops.h>

namespace lawa
{

template <typename Optype, typename Prec,
          typename T, typename Basis, typename Rhs>
unsigned
rank1AdaptiveAls(      Sepop<Optype>&                                A,
                       Prec&                                         P,
                       SepCoefficients<Lexicographical, T, Index1D>& v,
                       HTCoefficients<T, Basis>&                     x,
                       Rhs&                                          f,
                       std::vector<IndexSet<Index1D> >&              Lambda,
                 const Rank1AdaptiveAlsParams&                       params)
{
    // Basic input check
    assert(A.dim()  == (unsigned) x.dim());
    assert(A.dim()  == f.dim());
    assert(A.dim()  == Lambda.size());
    assert(v.rank() == 1);

    // ALS sweeps
    std::vector<IndexSet<Index1D> > start = Lambda;
    for (unsigned sweep=1; sweep<=params.max_sweep; ++sweep) {
        // Save result for stagnation check
        HTCoefficients<T, Basis> xold = x;

        for (unsigned j=1; j<=A.dim(); ++j) {
            // Compute projection (should be referenced to x)
            f.computeProjection(j);

            // Solve on leaf
            auto     v0    = v(1, j);
            Lambda[j-1]    = start[j-1];
            insert(x, v0, Lambda[j-1], j);
            unsigned numit = adaptiveLeaf(A, P, x, v0, Lambda, f, j,
                                          params.adaptiveLeaf);

            if (params.verbose) {
                std::cout << "rank1AdaptiveAls: Sweep " << sweep
                          << ", leaf " << j
                          << ", adaptiveLeaf required "
                          << numit << " iterations to reach r <= "
                          << params.adaptiveLeaf.tol << std::endl;
            }

            // Update x
            insert(x, v0, Lambda[j-1], j);
        }
        // Check stagnation
        xold.tree() = xold.tree() - x.tree();
        T stag      = nrm2(xold)/x.tree().L2norm();
        if (params.verbose) {
            std::cout << "rank1AdaptiveAls: Sweep " << sweep
                      << ", stagnation = " << stag << std::endl;
        }
        if (stag<=params.stag) {
            return sweep;
        }
    }

    std::cerr << "rank1AdaptiveAls: Reached max sweeps "
              << params.max_sweep << std::endl;
    return params.max_sweep;
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_RANK1ADAPTIVEALS_TCC
