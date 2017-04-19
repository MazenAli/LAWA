#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_CG_TCC
#define LAWA_METHODS_ADAPTIVE_SOLVERS_CG_TCC 1

#include <iostream>
#include <cmath>
#include <cassert>
#include <lawa/methods/adaptive/algorithms/coeffops.h>

namespace lawa
{

template <typename Optype, typename T, typename Basis>
unsigned
cg(      Sepop<Optype>&               A,
         flens::GeMatrix
         <flens::FullStorage
         <T, cxxblas::ColMajor> >&    U,
   const flens::GeMatrix
         <flens::FullStorage
         <T, cxxblas::ColMajor> >&    B,
   const flens::GeMatrix
         <flens::FullStorage
         <T, cxxblas::ColMajor> >&    Pj,
         HTCoefficients<T, Basis>&    u,
   const unsigned                     j,
   const IndexSet<Index1D>&           Lambda,
         T&                           res_cg,
   const T                            tol_cg,
   const unsigned                     maxit_cg,
   const bool                         verbose)
{
    assert(j>=1 && j<=A.dim());
    using flens::_;

    typedef typename flens::GeMatrix
                    <flens::FullStorage<T, cxxblas::ColMajor> >     Matrix;
    #ifdef CG_COMP_EVS
        typedef typename flens::DenseVector<flens::Array<T> >           Vector;
        typedef FLENS_DEFAULT_INDEXTYPE                                 INT;
        typedef typename flens::DenseVector<flens::Array<INT> >         iVector;

        /* Compute evs  */
        auto numr       = U.numRows();
        auto nume       = U.numCols();
        auto numc       = numr;
        Matrix Afull(numr, numc);
        for (int k=1; k<=numc; ++k) {
            Matrix ek(numr, nume);
            ek(k, 1)    = 1.;
            Matrix Aek;
            Aek         = redeval(A, ek, Pj, u, j, Lambda, Lambda);
            Afull(_, k) = Aek(_, 1);
        }

        Vector wr, wi, scale, rCondE, rCondV;
        Matrix vl, vr;
        INT ilo, ihi;
        T abNorm;
        Vector work;
        iVector iwork;
        flens::lapack::evx(flens::lapack::BALANCE::None,
                           false, false,
                           flens::lapack::SENSE::None,
                           Afull, wr, wi, vl, vr, ilo, ihi, scale, abNorm, rCondE,
                           rCondV, work, iwork);
        auto pos       = flens::blas::iamax(wi);
        T complex_part = wi(pos);
        T max          = wr(1);
        T min          = max;
        T absmax       = fabs(wr(1));
        T absmin       = absmax;
        for (int k=1; k<= wr.length(); ++k) {
            auto val    = wr(k);
            auto absval = fabs(val);
            max         = (val>max) ? val : max;
            absmax      = (absval>absmax) ? absval : absmax;
            min         = (val<min) ? val : min;
            absmin      = (absval<absmin) ? absval : absmin;
        }
        T cond         = absmax/absmin;
        std::cout << "cg: Complex part maxabs : " << complex_part << std::endl;
        std::cout << "cg: Max real part       : " << max          << std::endl;
        std::cout << "cg: Min real part       : " << min          << std::endl;
        std::cout << "cg: Condition           : " << cond         << std::endl;
        std::cout << "cg: epsilon             : " << Pj(1, 2)/Pj(1, 1) << "\n";
    #endif

    /* Initial residual */
    Matrix tmp = B-redeval(A, U, Pj, u, j, Lambda, Lambda);
    Matrix rk(tmp.numRows(), tmp.numCols());
    for (auto& lambda : Lambda) {
        auto i   = u.map()(lambda, j);
        rk(i, _) = tmp(i, _);
    }
    Matrix pk = rk, Apk(rk.numRows(), rk.numCols());
    T nrmb = flens::blas::nrm2(B.vectorView());
    T rho;
    T rho_old = flens::blas::dot(rk.vectorView(), rk.vectorView());
    res_cg         = std::sqrt(rho_old)/nrmb;
    if (verbose) {
        std::cout << "cg: Iteration " << 0 << " r = " << res_cg
                  << std::endl;
    }
    if (res_cg<=tol_cg) return 0;

    for (unsigned i=1; i<=maxit_cg; ++i) {
        /* Compute Apk */
        Apk  = redeval(A, pk, Pj, u, j, Lambda, Lambda);

        /* Compute corrections */
        T pAp = flens::blas::dot(pk.vectorView(), Apk.vectorView());
        T ak  = rho_old/pAp;
        U    += ak*pk;
        rk   -= ak*Apk;
        rho   = flens::blas::dot(rk.vectorView(), rk.vectorView());

        res_cg     = std::sqrt(rho);
        res_cg    /= nrmb;
        if (verbose) {
            std::cout << "cg: Iteration " << i << " residual = " << res_cg
                      << std::endl;
        }
        if (res_cg<=tol_cg) return i;

        /* Update */
        T bk       = rho/rho_old;
        rho_old    = rho;
        pk        *= bk;
        pk        += rk;
    }

    std::cerr << "cg: Reached max iterations " << maxit_cg << "\n";
    return maxit_cg;
}

template <typename Optype, typename T, typename Basis>
unsigned
cg(      Sepop<Optype>&               A,
         Sepdiagscal<Basis>&          S,
         flens::GeMatrix
         <flens::FullStorage
         <T, cxxblas::ColMajor> >&    U,
   const flens::GeMatrix
         <flens::FullStorage
         <T, cxxblas::ColMajor> >&    B,
   const flens::GeMatrix
         <flens::FullStorage
         <T, cxxblas::ColMajor> >&    Pj,
         HTCoefficients<T, Basis>&    u,
   const unsigned                     j,
   const IndexSet<Index1D>&           Lambda,
         T&                           res_cg,
   const T                            tol_cg,
   const unsigned                     maxit_cg)
{
    assert(j>=1 && j<=A.dim());

    typedef typename flens::GeMatrix
                    <flens::FullStorage<T, cxxblas::ColMajor> >     Matrix;

    #ifdef CG_COMP_EVS
        typedef typename flens::DenseVector<flens::Array<T> >           Vector;
        typedef FLENS_DEFAULT_INDEXTYPE                                 INT;
        typedef typename flens::DenseVector<flens::Array<INT> >         iVector;

        /* Compute evs  */
        using flens::_;
        auto numr       = U.numRows();
        auto nume       = U.numCols();
        auto numc       = numr;
        Matrix Afull(numr, numc);
        for (int k=1; k<=numc; ++k) {
            Matrix ek(numr, nume);
            ek(k, 1)    = 1.;
            Matrix Aek;
            Aek         = redeval(A, P, ek, Pj, u, j, Lambda, Lambda);
            Afull(_, k) = Aek(_, 1);
        }

        Vector wr, wi, scale, rCondE, rCondV;
        Matrix vl, vr;
        INT ilo, ihi;
        T abNorm;
        Vector work;
        iVector iwork;
        flens::lapack::evx(flens::lapack::BALANCE::None,
                           false, false,
                           flens::lapack::SENSE::None,
                           Afull, wr, wi, vl, vr, ilo, ihi, scale, abNorm, rCondE,
                           rCondV, work, iwork);
        auto pos       = flens::blas::iamax(wi);
        T complex_part = wi(pos);
        T max          = wr(1);
        T min          = max;
        T absmax       = fabs(wr(1));
        T absmin       = absmax;
        for (int k=1; k<= wr.length(); ++k) {
            auto val    = wr(k);
            auto absval = fabs(val);
            max         = (val>max) ? val : max;
            absmax      = (absval>absmax) ? absval : absmax;
            min         = (val<min) ? val : min;
            absmin      = (absval<absmin) ? absval : absmin;
        }
        T cond         = absmax/absmin;
        std::cout << "cg: Complex part maxabs : " << complex_part << std::endl;
        std::cout << "cg: Max real part       : " << max          << std::endl;
        std::cout << "cg: Min real part       : " << min          << std::endl;
        std::cout << "cg: Condition           : " << cond         << std::endl;
    #endif

    /* Initial residual */
    Matrix rk = B-redeval(A, S, U, Pj, u, j, Lambda, Lambda);
    Matrix pk = rk, Apk;
    T nrmb = flens::blas::nrm2(B.vectorView());
    T rho;
    T rho_old = flens::blas::dot(rk.vectorView(), rk.vectorView());
    res_cg         = std::sqrt(rho_old)/nrmb;
    #ifdef VERBOSE
        std::cout << "cg: Iteration " << 0 << " r = " << res_cg
                  << std::endl;
    #endif
    if (res_cg<=tol_cg) return 0;

    for (unsigned i=1; i<=maxit_cg; ++i) {
        /* Compute Apk */
        Apk  = redeval(A, S, pk, Pj, u, j, Lambda, Lambda);

        /* Compute corrections */
        T pAp = flens::blas::dot(pk.vectorView(), Apk.vectorView());
        T ak  = rho_old/pAp;
        U    += ak*pk;
        rk   -= ak*Apk;
        rho   = flens::blas::dot(rk.vectorView(), rk.vectorView());

        res_cg     = std::sqrt(rho);
        res_cg    /= nrmb;
        #ifdef VERBOSE
            std::cout << "cg: Iteration " << i << " residual = " << res_cg
                      << std::endl;
        #endif
        if (res_cg<=tol_cg) return i;

        /* Update */
        T bk       = rho/rho_old;
        rho_old    = rho;
        pk        *= bk;
        pk        += rk;
    }

    std::cerr << "cg: Reached max iterations " << maxit_cg << "\n";
    return maxit_cg;
}


template <typename Optype, typename Prec, typename T, typename Basis>
unsigned
cg(      Sepop<Optype>&               A,
         Prec&                        P,
         flens::GeMatrix
         <flens::FullStorage
         <T, cxxblas::ColMajor> >&    U,
   const flens::GeMatrix
         <flens::FullStorage
         <T, cxxblas::ColMajor> >&    B,
   const flens::GeMatrix
         <flens::FullStorage
         <T, cxxblas::ColMajor> >&    Pj,
         HTCoefficients<T, Basis>&    u,
   const unsigned                     j,
   const IndexSet<Index1D>&           Lambda,
         T&                           res_cg,
   const T                            tol_cg,
   const unsigned                     maxit_cg)
{
    assert(j>=1 && j<=A.dim());

    typedef typename flens::GeMatrix
                    <flens::FullStorage<T, cxxblas::ColMajor> >     Matrix;

    #ifdef CG_COMP_EVS
        typedef typename flens::DenseVector<flens::Array<T> >           Vector;
        typedef FLENS_DEFAULT_INDEXTYPE                                 INT;
        typedef typename flens::DenseVector<flens::Array<INT> >         iVector;

        /* Compute evs  */
        using flens::_;
        auto numr       = U.numRows();
        auto nume       = U.numCols();
        auto numc       = numr;
        Matrix Afull(numr, numc);
        for (int k=1; k<=numc; ++k) {
            Matrix ek(numr, nume);
            ek(k, 1)    = 1.;
            Matrix Aek;
            Aek         = redeval(A, P, ek, Pj, u, j, Lambda, Lambda);
            Afull(_, k) = Aek(_, 1);
        }

        Vector wr, wi, scale, rCondE, rCondV;
        Matrix vl, vr;
        INT ilo, ihi;
        T abNorm;
        Vector work;
        iVector iwork;
        flens::lapack::evx(flens::lapack::BALANCE::None,
                           false, false,
                           flens::lapack::SENSE::None,
                           Afull, wr, wi, vl, vr, ilo, ihi, scale, abNorm, rCondE,
                           rCondV, work, iwork);
        auto pos       = flens::blas::iamax(wi);
        T complex_part = wi(pos);
        T max          = wr(1);
        T min          = max;
        T absmax       = fabs(wr(1));
        T absmin       = absmax;
        for (int k=1; k<= wr.length(); ++k) {
            auto val    = wr(k);
            auto absval = fabs(val);
            max         = (val>max) ? val : max;
            absmax      = (absval>absmax) ? absval : absmax;
            min         = (val<min) ? val : min;
            absmin      = (absval<absmin) ? absval : absmin;
        }
        T cond         = absmax/absmin;
        std::cout << "cg: Complex part maxabs : " << complex_part << std::endl;
        std::cout << "cg: Max real part       : " << max          << std::endl;
        std::cout << "cg: Min real part       : " << min          << std::endl;
        std::cout << "cg: Condition           : " << cond         << std::endl;
        std::cout << "cg: epsilon             : " << Pj(1, 2)/Pj(1, 1) << "\n";
    #endif

    /* Adjust for preconditioning */
    Matrix DB = prec(P, B, u, j, Lambda);
    U         = remove_prec(P, U, u, j, Lambda);

    /* Initial residual */
    Matrix rk = DB-redeval(A, P, U, Pj, u, j, Lambda, Lambda);
    Matrix pk = rk, Apk;
    T nrmb    = flens::blas::nrm2(DB.vectorView());
    T rho;
    T rho_old = flens::blas::dot(rk.vectorView(), rk.vectorView());
    res_cg    = std::sqrt(rho_old)/nrmb;
    #ifdef VERBOSE
        std::cout << "cg: Iteration " << 0 << " r = " << res_cg
                  << std::endl;
    #endif
    if (res_cg<=tol_cg) return 0;

    for (unsigned i=1; i<=maxit_cg; ++i) {
        /* Compute Apk */
        Apk  = redeval(A, P, pk, Pj, u, j, Lambda, Lambda);

        /* Compute corrections */
        T pAp = flens::blas::dot(pk.vectorView(), Apk.vectorView());
        T ak  = rho_old/pAp;
        U    += ak*pk;
        rk   -= ak*Apk;
        rho   = flens::blas::dot(rk.vectorView(), rk.vectorView());

        res_cg = std::sqrt(rho)/nrmb;
        #ifdef VERBOSE
            std::cout << "cg: Iteration " << i << " residual = " << res_cg
                      << std::endl;
        #endif
        if (res_cg<=tol_cg) {
            U = prec(P, U, u, j, Lambda);
            return i;
        }

        /* Update */
        T bk       = rho/rho_old;
        rho_old    = rho;
        pk        *= bk;
        pk        += rk;
    }

    U = prec(P, U, u, j, Lambda);
    std::cerr << "cg: Reached max iterations " << maxit_cg << "\n";
    return maxit_cg;
}


template <typename Optype, typename Prec, typename T, typename Basis>
unsigned
cg_rank1prec(    Sepop<Optype>&               A,
                 Prec&                        P,
                 flens::GeMatrix
                 <flens::FullStorage
                 <T, cxxblas::ColMajor> >&    U,
           const flens::GeMatrix
                 <flens::FullStorage
                 <T, cxxblas::ColMajor> >&    B,
           const flens::GeMatrix
                 <flens::FullStorage
                 <T, cxxblas::ColMajor> >&    Pj,
                 HTCoefficients<T, Basis>&    u,
           const unsigned                     j,
           const IndexSet<Index1D>&           Lambda,
                 T&                           res_cg,
           const T                            tol_cg,
           const unsigned                     maxit_cg,
           const bool                         verbose)
{
    assert(j>=1 && j<=A.dim());
    using flens::_;

    typedef typename flens::GeMatrix
                    <flens::FullStorage<T, cxxblas::ColMajor> >     Matrix;

    #ifdef CG_COMP_EVS
        typedef typename flens::DenseVector<flens::Array<T> >           Vector;
        typedef FLENS_DEFAULT_INDEXTYPE                                 INT;
        typedef typename flens::DenseVector<flens::Array<INT> >         iVector;

        /* Compute evs  */
        auto numr = U.numRows();
        Matrix Afull(numr, numr);
        for (unsigned k=1; k<=numr; ++k) {
            Matrix ek(numr, U.numCols());
            ek(k, 1)      = 1.;
            Matrix Aek    = redeval(A, ek, Pj, u, j, Lambda, Lambda);
            Aek           = apply_precsq(P, Aek);
            Afull(_, k)   = Aek(_, 1);
        }
        Vector wr, wi, scale, rCondE, rCondV;
        Matrix vl, vr;
        INT ilo, ihi;
        T abNorm;
        Vector work;
        iVector iwork;
        flens::lapack::evx(flens::lapack::BALANCE::None,
                           false, false,
                           flens::lapack::SENSE::None,
                           Afull, wr, wi, vl, vr, ilo, ihi, scale, abNorm, rCondE,
                           rCondV, work, iwork);
        auto pos       = flens::blas::iamax(wi);
        T complex_part = wi(pos);
        T max          = wr(1);
        T min          = max;
        T absmax       = fabs(wr(1));
        T absmin       = absmax;
        for (int k=1; k<= wr.length(); ++k) {
            auto val    = wr(k);
            auto absval = fabs(val);
            max         = (val>max) ? val : max;
            absmax      = (absval>absmax) ? absval : absmax;
            min         = (val<min) ? val : min;
            absmin      = (absval<absmin) ? absval : absmin;
        }
        T cond         = absmax/absmin;
        std::cout << "cg: Complex part maxabs : " << complex_part << std::endl;
        std::cout << "cg: Max real part       : " << max          << std::endl;
        std::cout << "cg: Min real part       : " << min          << std::endl;
        std::cout << "cg: Condition           : " << cond         << std::endl;
        std::cout << "cg: epsilon             : " << Pj(1, 2)/Pj(1, 1) << "\n";
    #endif

    /* Initial residual */
    Matrix tmp = B-redeval(A, U, Pj, u, j, Lambda, Lambda);
    Matrix rk(tmp.numRows(), tmp.numCols());
    for (auto& lambda : Lambda) {
        auto i   = u.map()(lambda, j);
        rk(i, _) = tmp(i, _);
    }

    Matrix zk = apply_precsq(P, rk);
    Matrix pk = zk, Apk(rk.numRows(), rk.numCols());
    T nrmb    = flens::blas::nrm2(B.vectorView());
    T rho;
    T rho_old = flens::blas::dot(rk.vectorView(), zk.vectorView());
    res_cg    = flens::blas::nrm2(rk.vectorView())/nrmb;
    if (verbose) {
        std::cout << "cg: Iteration " << 0 << " r = " << res_cg
                  << std::endl;
    }
    if (res_cg<=tol_cg) return 0;

    for (unsigned i=1; i<=maxit_cg; ++i) {
        /* Compute Apk */
        Apk  = redeval(A, pk, Pj, u, j, Lambda, Lambda);

        /* Compute corrections */
        T pAp = flens::blas::dot(pk.vectorView(), Apk.vectorView());
        T ak  = rho_old/pAp;
        U    += ak*pk;
        rk   -= ak*Apk;
        zk    = apply_precsq(P, rk);
        rho   = flens::blas::dot(rk.vectorView(), zk.vectorView());

        res_cg = flens::blas::nrm2(rk.vectorView())/nrmb;
        if (verbose) {
            std::cout << "cg: Iteration " << i << " residual = " << res_cg
                      << std::endl;
        }
        if (res_cg<=tol_cg) return i;

        /* Update */
        T bk     = rho/rho_old;
        rho_old  = rho;
        pk      *= bk;
        pk      += zk;
    }

    std::cerr << "cg: Reached max iterations " << maxit_cg << "\n";
    return maxit_cg;
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_CG_TCC
