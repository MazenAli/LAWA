#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_SPLITTING_TCC
#define LAWA_METHODS_ADAPTIVE_SOLVERS_SPLITTING_TCC 1

namespace lawa
{

template <typename Optype, typename Prec, typename T, typename Basis>
unsigned
symm_splitting(Sepop<Optype>&             A,
             Prec&                        P,
             flens::GeMatrix
             <flens::FullStorage
             <T, cxxblas::ColMajor> >&    V,
       const flens::GeMatrix
             <flens::FullStorage
             <T, cxxblas::ColMajor> >&    B,
       const flens::GeMatrix
             <flens::FullStorage
             <T, cxxblas::ColMajor> >&    Pj,
             HTCoefficients<T, Basis>&    u,
       const unsigned                     j,
       const IndexSet<Index1D>&           Lambda,
             T&                           res,
       const T                            tol,
       const unsigned                     maxit)
{
    assert(j>=1 && j<=A.dim());

    using flens::_;
    typedef typename flens::GeMatrix
                    <flens::FullStorage<T, cxxblas::ColMajor> >     Matrix;

    Matrix VG1, VG2, DB, AVG2, rk, SS;
    const Matrix& G1 = Pj(_, _(1, V.numCols()));
    const Matrix& G2 = Pj(_, _(V.numCols()+1, Pj.numCols()));

    /* Adjust for scaling */
    V = remove_prec(P, V, u, j, Lambda);

    /* Initial residual */
    DB   = prec(P, B, u, j, Lambda);
    flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans, 1., V, G1, 0., VG1);
    flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans, 1., V, G2, 0., VG2);
    VG1  = precsq(P, VG1, u, j, Lambda);
    AVG2 = eval(A, P, VG2, u, j, Lambda, Lambda);
    rk   = DB-AVG2-VG1;

    res     = flens::blas::nrm2(rk.vectorView());
    T nrmb  = flens::blas::nrm2(DB.vectorView());
    res    /= nrmb;
    #ifdef VERBOSE
        std::cout << "symm_splitting: Iteration " << 0 << " r = " << res
                  << std::endl;
    #endif
    if (res<=tol) {
        V = prec(P, V, u, j, Lambda);
        return 0;
    }

    SS = DB;
    flens::blas::scal(0.5, SS);
    Matrix c;
    for (unsigned i=1; i<=maxit; ++i) {
        /* First stage */
        c  = prec(P, V, u, j, Lambda);
        c -= redeval(A, P, V, Pj, u, j, Lambda, Lambda);
        flens::blas::scal(-0.5, c);
        c += SS;
        T res_cg, tol_cg = 1e-08, maxit_cg = 100;
        unsigned it = cg(A, P, V, B, Pj, u, j,
                         Lambda,
                         res_cg, tol_cg, maxit_cg);
        #ifdef VERBOSE
            std::cout << "symm_splitting: 1st stage, CG required " << it
                      << " iterations to reach "
                      << res_cg << "\n";
        #endif

        /* Second stage */
        c  = prec(P, V, u, j, Lambda);
        c -= redeval(A, P, V, Pj, u, j, Lambda, Lambda);
        flens::blas::scal(-0.5, c);
        c += SS;
        c  = remove_precsq(P, c, u, j, Lambda);
        V  = c;
        flens::blas::scal(1./Pj(1, 1), V);

        /* Third stage */
        c  = prec(P, V, u, j, Lambda);
        c -= redeval(A, P, V, Pj, u, j, Lambda, Lambda);
        flens::blas::scal(-0.5, c);
        c += SS;
        it = cg(A, P, V, B, Pj, u, j,
                Lambda,
                res_cg, tol_cg, maxit_cg);
        #ifdef VERBOSE
            std::cout << "symm_splitting: 3rd stage, CG required " << it
                      << " iterations to reach "
                      << res_cg << "\n";
        #endif

        flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans, 1., V, G1, 0., VG1);
        flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans, 1., V, G2, 0., VG2);
        VG1  = precsq(P, VG1, u, j, Lambda);
        AVG2 = eval(A, P, VG2, u, j, Lambda, Lambda);
        rk   = DB-AVG2-VG1;

        res     = flens::blas::nrm2(rk.vectorView());
        res    /= nrmb;
        #ifdef VERBOSE
            std::cout << "symm_splitting: Iteration " << i << " r = " << res
                      << std::endl;
        #endif
        if (res<=tol) {
            V = prec(P, V, u, j, Lambda);
            return i;
        }
    }

    /* Adjust for scaling */
    V = prec(P, V, u, j, Lambda);
    std::cerr << "symm_splitting: Reached max iterations " << maxit << "\n";
    return maxit;
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_SPLITTING_TCC
