#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_WAVPCG_TCC
#define LAWA_METHODS_ADAPTIVE_SOLVERS_WAVPCG_TCC 1

#include <iostream>
#include <cassert>
#include <cmath>

namespace lawa
{

template <typename Operator, typename Preconditioner, typename T>
unsigned
wavPcg(      Operator&                                  A,
             Preconditioner&                            P,
             Coefficients<Lexicographical, T, Index1D>& x,
       const Coefficients<Lexicographical, T, Index1D>& b,
             T&                                         res_cg,
       const T                                          tol,
       const unsigned                                   maxit,
       const bool                                       verbose)
{
    assert(tol>0);

    Coefficients<Lexicographical, T, Index1D> r;
    Coefficients<Lexicographical, T, Index1D> res;
    Coefficients<Lexicographical, T, Index1D> z;
    Coefficients<Lexicographical, T, Index1D> p;
    Coefficients<Lexicographical, T, Index1D> Ap;

    // Initial residual
    r   = b - A(x);
    res = r;
    P(res);
    z   = res;
    P(z);
    p = z;
    T rho;
    T rho_old = z*r;
    res_cg    = res.norm(2.);

    if (verbose) {
        std::cout << "wavPcg: Iteration 0, r = " << res_cg << std::endl;
    }
    if (res_cg<=tol) return 0;

    // pcg iterations
    for (unsigned k=1; k<=maxit; ++k) {
        // Compute Ap
        Ap = A(p);

        // Compute corrections
        T pAp = p*Ap;
        T a   = rho_old/pAp;
        x    += a*p;
        r    -= a*Ap;

        // Check residual
        res    = r;
        P(res);
        z      = res;
        P(z);
        rho    = z*r;
        res_cg = res.norm(2.);
        if (verbose) {
            std::cout << "wavPcg: Iteration " << k << ", r = " << res_cg
                      << std::endl;
        }
        if (res_cg<=tol) return k;

        // Update
        T b     = rho/rho_old;
        rho_old = rho;
        p       = z+b*p;
    }

    std::cerr << "wavPcg: Reached max iterations " << maxit << std::endl;
    return maxit;
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_WAVPCG_TCC
