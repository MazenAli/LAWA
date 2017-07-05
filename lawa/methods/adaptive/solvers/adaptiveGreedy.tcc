#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_ADAPTIVEGREEDY_TCC
#define LAWA_METHODS_ADAPTIVE_SOLVERS_ADAPTIVEGREEDY_TCC 1

#include <cassert>

namespace lawa
{

unsigned
adaptiveGreedy(Engine* ep,
                Sepop<Optype>&  A,
                Sepdiagscal<Basis>& S,
                Prec&               P,
                HTCoefficients<T, Basis> u)
{
    // Basic input check
    assert(A.dim() == S.dim());
    assert(A.dim() == (unsigned) u.dim());
    assert(A.dim() == f.dim());
    assert(A.dim() == Lambda.size());

    // Typedefs
    typedef flens::GeMatrix<flens::FullStorage<T, flens::ColMajor> >   Matrix;
    typedef flens::DenseVector<flens::Array<FLENS_DEFAULT_INDEXTYPE> > IVector;
    typedef std::vector<Matrix>                                        MatVec;

    // Aux variables
    HTCoefficients<T, Basis>                     r(u.dim(), u.basis(), u.map());
    HTCoefficients<T, Basis>                     F(u.dim(), u.basis(), u.map());
    r.tree().set_tree(u.tree());
    F.tree().set_tree(u.tree());
    SepCoefficients<Lexicographical, T, Index1D> Fcp(f.rank(), f.dim());

    std::vector<IndexSet<Index1D> >              total   = Lambda;
    std::vector<IndexSet<Index1D> >              current = Lambda;

    T trunc = 1e-08;

    // Initial RHS
    genCoefficients(Fcp, f, Lambda);
    set(F, Fcp, Lambda);

    // Initial residual
    r        = eval(A, u, Lambda, Lambda);
    r.tree() = F.tree()-r.tree();
    r        = applyScaleTT(S, F, Lambda, trunc);
    T res    = nrm2(r);

    if (params.verbose) {
        std::cout << "adaptiveGreedy: Update 0, r = " << res << std::endl;
    }
    if (res<=params.tol || params.maxit<1) return 0;

    // Rank 1 ALS
    SepCoefficients<Lexicographical, T, Index1D> v(1, A.dim());
    genCoefficientsRnd(v, Lambda, 1e-03, 1);
    ProjRhs<T, Basis> fj(Fcp, current, f, u);
    unsigned numit = rank1AdaptiveAls(A, P, v, u, fj, current, params.r1Als);
    Lambda  = unify(Lambda, current);

    // Residual
    r   = szoneres(A, u, F, Fcp, f, Lambda, Lambda, total);
    r   = applyScaleTT(S, F, total, trunc);
    res = nrm2(r);

    if (params.verbose) {
        std::cout << "adaptiveGreedy: Update 1, adaptiveAls required "
                  << numit << " sweeps, r = " << res << std::endl;
    }
    if (res<=params.tol || params.maxit<2) return 1;

    // Updates
    for (unsigned k=2; k<=params.maxit) {
        // Rank 1 update
        genCoefficientsRnd(v, Lambda, 1e-03, 1);
        HTCoefficients<T, Basis> x(u.dim(), u.basis(), u.map());
        x.tree().set_tree(u.tree());
        set(x, v, Lambda);
        ProjResidual<T, Basis, Optype f_Auj(Fcp, f, x, u, A, current, Lambda);
        current  = Lambda;
        numit    = rank1AdaptiveAls(A, P, v, u, f_Auj, current, params.r1Als);
        u.tree() = x.tree()+u.tree();
        Lambda   = unify(Lambda, current);

        // DMRG for Galerkin solve on reduced space
        u.tree().orthogonalize();
        genCoefficients(Fcp, f, Lambda);
        set(F, Fcp, Lambda);
        MatVec fred = reduce_rhs(u, F);
        MatVec Ared = reduce(A, u, Lambda, Lambda);
        std::vector<Matrix> x0(fred.size()-1);
        IVector ranks(x0.size());
        htucker::extract_core(u.tree(), x0, ranks);
        x0 = optTTcoreLaplace(ep, Ared, fred, x0, ranks, params.dmrgTTcore);
        htucker::insert_core(x.tree(), x0, ranks);

        // Residual
        r   = szoneres(A, u, F, Fcp, f, Lambda, Lambda, total);
        r   = applyScaleTT(S, F, total, trunc);
        res = nrm2(r);

        // Check
        if (params.verbose) {
            std::cout << "adaptiveGreedy: Update " << k
                      << ",  adaptiveAls required "
                      << numit << " sweeps, r = " << res << std::endl;
        }
        if (res<=params.tol) return k;
    }

    std::cerr << "adaptiveGreedy: Reached max updates " << params.maxit
              << std::endl;
    return params.maxit;
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_ADAPTIVEGREEDY_TCC
