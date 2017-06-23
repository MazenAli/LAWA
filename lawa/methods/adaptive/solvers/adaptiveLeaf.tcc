#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_adaptiveLeaf_TCC
#define LAWA_METHODS_ADAPTIVE_SOLVERS_adaptiveLeaf_TCC 1

#include <cassert>

#include <htucker/htucker.h>

#include <lawa/methods/adaptive/algorithms/coeffops.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>

namespace lawa
{

template <typename Optype, typename T, typename Prec, typename Rhs>
unsigned
adaptiveLeaf_laplace(      Sepop<Optype>&                              A,
                     const flens::GeMatrix<
                           flens::FullStorage< T, flens::ColMajor > >& PAj,
                           Prec&                                       P,
                           HTCoefficients<T, Basis>&                   u,
                           IndexSet<Index1D>&                          Lambda,
                     const SeparableRHSD<T, Basis>&                    f,
                     const HTCoefficients<T, Basis>&                   b,
                     const flens::GeMatrix<
                           flens::FullStorage< T, flens::ColMajor > >& Pbj,
                     const unsigned                                    j,
                     const AdaptiveAlsParams&                          p);
{
    // Basic input check
    assert(A.dim() == (unsigned) u.dim());
    assert(A.dim() == (unsigned) f.dim());
    assert(A.dim() == (unsigned) b.dim());
    assert(Lambda.size()>0);
    assert(j>=1 && j<=A.dim());

    typedef typename flens::GeMatrix<
                     flens::FullStorage< T, flens::ColMajor > >     Matrix;

    // Extract factor
    Coefficients<Lexicographical, T, Index> v = extract(u, Lambda, j);

    // Evaluate rhs
    // 1. Extract leaf
    // 2. Contract
    // 3. Convert to coefficients
    // 4. Galerkin solve
    // 5. Extend to new index set
    // 6. Evaluate rank_j(f) rhs's on the extended difference
    // 7. Convert to integer indices and add to leaf
    // 8. Contract. Proceed with Galerkin solve etc.
    // 9. Insert final result into the rhs tensor
    //
    // Implement separately cg Galerkin solve
    // Implement redeval in simple form here
    htucker::DimensionIndex idx(1);
    idx[0]    = j;
    Matrix UB = extract(b.tree(), idx);
    Matrix B;
    flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans, 1.0, UB, Pbj, 0.0, B);
    Coefficients<Lexicographical, T, Index> r;
    convert(B, r, u);

    // Initial residual
    r = r-redeval(A, v)

    // AWGM iterations
    for ()
        // Galerkin solve
        cg <- new routine here

        // Extend index set

        // Evaluate residual

        // Stopping criterion
            // Insert factor

        // Bulk chasing

    return;
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_adaptiveLeaf_TCC
