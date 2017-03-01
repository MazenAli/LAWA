#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_LEAF_SOLVE_TCC
#define LAWA_METHODS_ADAPTIVE_SOLVERS_LEAF_SOLVE_TCC 1

namespace lawa
{

template <typename Optype, typename Prec, typename T, typename Basis>
unsigned
leaf_solve(Sepop<Optype>&               A,
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

    Matrix VG1, VG2, DB, AVG2, rk;
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
        std::cout << "leaf_solve: Iteration " << 0 << " r = " << res
                  << std::endl;
    #endif
    if (res<=tol) {
        V = prec(P, V, u, j, Lambda);
        return 0;
    }

    for (unsigned i=1; i<=maxit; ++i) {
        /* Stage 1 */
        Matrix defect = DB-AVG2;
        defect = remove_precsq(P, defect, u, j, Lambda);
        // Solve
        flens::blas::scal(1./G1(1,1), defect);

        V      = defect;

        flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans, 1., V, G1, 0., VG1);
        flens::blas::mm(cxxblas::NoTrans, cxxblas::Trans, 1., V, G2, 0., VG2);
        VG1  = precsq(P, VG1, u, j, Lambda);
        AVG2 = eval(A, P, VG2, u, j, Lambda, Lambda);
        rk   = DB-AVG2-VG1;

        res_fp  = flens::blas::nrm2(rk.vectorView());
        T nrmb  = flens::blas::nrm2(DB.vectorView());
        res_fp /= nrmb;
        #ifdef VERBOSE
            std::cout << "leaf_solve: Iteration " << i << " r = " << res_fp
                      << std::endl;
        #endif
        if (res<=tol) {
            V = prec(P, V, u, j, Lambda);
            return i;
        }
    }

    /* Adjust for scaling */
    V = prec(P, V, u, j, Lambda);
    std::cerr << "leaf_solve: Reached max iterations " << maxit << "\n";
    return maxit;
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_LEAF_SOLVE_TCC
