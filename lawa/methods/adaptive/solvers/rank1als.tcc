#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_RANK1ALS_TCC
#define LAWA_METHODS_ADAPTIVE_SOLVERS_RANK1ALS_TCC 1

#include <iostream>
#include <cmath>
#include <cassert>
#include <htucker/htucker.h>
#include <flens/flens.cxx>
#include <lawa/methods/adaptive/solvers/cg.h>
#include <lawa/methods/adaptive/solvers/splitting.h>
#include <lawa/preconditioners/preconditioners1d/H1normpreconditioner1d.h>
#include <lawa/methods/adaptive/algorithms/coeffops.h>

namespace lawa
{

template <typename Optype, typename T, typename Basis>
unsigned
rank1als_sym(       Sepop<Optype>&                      A,
                    HTCoefficients<T, Basis>&           x,
              const HTCoefficients<T, Basis>&           b,
              const std::vector<IndexSet<Index1D> >&    Lambda,
                    T&                                  residual,
              const bool                                orthog,
              const T                                   tol,
              const unsigned                            max_sweep,
              const T                                   tol_cg,
              const unsigned                            maxit_cg)
{
    typedef typename flens::GeMatrix
                    <flens::FullStorage<T, cxxblas::ColMajor> >     Matrix;
    typedef typename htucker::HTuckerTree<T>                        HTTree;
    typedef HTCoefficients<T, Basis>                                HTCoeff;

    /* Initial residual */
    HTCoeff Ax;
    T nrmb = nrm2(const_cast<HTCoeff&>(b));

    Ax         = eval(A, x, Lambda, Lambda);
    HTTree res = b.tree()-Ax.tree();
    res.orthogonalize();
    residual   = res.L2normorthogonal()/nrmb;
    #ifdef VERBOSE
        std::cout << "rank1als_sym: Sweep " << 0 << ", leaf " << 0
                  << ", residual = " << residual << "\n";
    #endif
    if (residual<=tol) return 0;

    /* ALS sweeps */
    for (unsigned sweep=1; sweep<=max_sweep; ++sweep) {
        for (unsigned j=1; j<=A.dim(); ++j) {
            /* Compute projection */
            Matrix Pj  = projection(Ax.tree(), x.tree(), j);
            Matrix B   = contract(b.tree(), x.tree(), j);

            /* CG solve on leaf*/
            htucker::DimensionIndex idx(1);
            idx[0]           = j;
            Matrix   xk      = extract(x.tree(), idx);
            T res_cg;
            unsigned it      = cg(A, xk, B, Pj, x, j,
                                  Lambda[j-1], res_cg, tol_cg, maxit_cg);
            #ifdef VERBOSE
                std::cout << "rank1als_sym: CG required " << it
                          << " iterations to reach "
                          << res_cg << "\n";
            #endif

            /* Update x */
            insert(x.tree(), xk, idx);
            if (orthog) x.tree().orthogonalize();

            /* Full residual */
            Ax         = eval(A, x, Lambda, Lambda);
            res        = b.tree()-Ax.tree();
            res.orthogonalize();
            residual   = res.L2normorthogonal()/nrmb;
            #ifdef VERBOSE
                std::cout << "rank1als_sym: Sweep " << sweep << ", leaf " << j
                          << ", residual = " << residual << "\n";
            #endif
            if (residual<=tol) return sweep;
        }
    }

    std::cerr << "rank1als_sym: Reached max sweeps " << max_sweep << "\n";
    return max_sweep;
}


template <typename Optype, typename T, typename Basis>
unsigned
rank1als_sym(       Sepop<Optype>&                      A,
                    Sepdiagscal<Basis>&                 S,
                    HTCoefficients<T, Basis>&           x,
              const HTCoefficients<T, Basis>&           b,
              const std::vector<IndexSet<Index1D> >&    Lambda,
                    T&                                  residual,
              const bool                                orthog,
              const T                                   tol,
              const unsigned                            max_sweep,
              const T                                   tol_cg,
              const unsigned                            maxit_cg)
{
    typedef typename flens::GeMatrix
                    <flens::FullStorage<T, cxxblas::ColMajor> >     Matrix;
    typedef typename htucker::HTuckerTree<T>                        HTTree;
    typedef HTCoefficients<T, Basis>                                HTCoeff;

    /* Initial residual */
    HTCoeff Ax, Sb = fixeval_notrunc(S, const_cast<HTCoeff&>(b), Lambda);
    T nrmb = nrm2(Sb);

    Ax         = fixeval_notrunc(S, x, Lambda);
    Ax         = eval(A, Ax, Lambda, Lambda);
    Ax         = fixeval_notrunc(S, Ax, Lambda);
    HTTree res = Sb.tree()-Ax.tree();
    res.orthogonalize();
    residual   = res.L2normorthogonal()/nrmb;
    #ifdef VERBOSE
        std::cout << "rank1als_sym: Sweep " << 0 << ", leaf " << 0
                  << ", residual = " << residual << "\n";
    #endif
    if (residual<=tol) return 0;

    /* ALS sweeps */
    for (unsigned sweep=1; sweep<=max_sweep; ++sweep) {
        for (unsigned j=1; j<=A.dim(); ++j) {
            /* Compute projection */
            Matrix Pj  = projection(Ax.tree(), x.tree(), j);
            Matrix B   = contract(Sb.tree(), x.tree(), j);

            /* CG solve on leaf*/
            htucker::DimensionIndex idx(1);
            idx[0]           = j;
            Matrix   xk      = extract(x.tree(), idx);
            T res_cg;
            unsigned it      = cg(A, S, xk, B, Pj, x, j,
                                  Lambda[j-1], res_cg, tol_cg, maxit_cg);
            #ifdef VERBOSE
                std::cout << "rank1als_sym: CG required " << it
                          << " iterations to reach "
                          << res_cg << "\n";
            #endif

            /* Update x */
            insert(x.tree(), xk, idx);
            if (orthog) x.tree().orthogonalize();

            /* Full residual */
            Ax         = fixeval_notrunc(S, x, Lambda);
            Ax         = eval(A, Ax, Lambda, Lambda);
            Ax         = fixeval_notrunc(S, Ax, Lambda);
            res        = Sb.tree()-Ax.tree();
            res.orthogonalize();
            residual   = res.L2normorthogonal()/nrmb;
            #ifdef VERBOSE
                std::cout << "rank1als_sym: Sweep " << sweep << ", leaf " << j
                          << ", residual = " << residual << "\n";
            #endif
            if (residual<=tol) return sweep;
        }
    }

    std::cerr << "rank1als_sym: Reached max sweeps " << max_sweep << "\n";
    return max_sweep;
}


template <typename Optype, typename Prec, typename T, typename Basis>
unsigned
rank1als_sym(       Sepop<Optype>&                      A,
                    Prec&                               P,
                    HTCoefficients<T, Basis>&           x,
              const HTCoefficients<T, Basis>&           b,
              const std::vector<IndexSet<Index1D> >&    Lambda,
                    T&                                  residual,
              const bool                                orthog,
              const T                                   tol,
              const unsigned                            max_sweep,
              const T                                   tol_cg,
              const unsigned                            maxit_cg)
{
    typedef typename flens::GeMatrix
                    <flens::FullStorage<T, cxxblas::ColMajor> >     Matrix;
    typedef typename htucker::HTuckerTree<T>                        HTTree;
    typedef HTCoefficients<T, Basis>                                HTCoeff;

    /* Initial residual */
    HTCoeff Ax;
    T nrmb = nrm2(const_cast<HTCoeff&>(b));

    Ax         = eval(A, x, Lambda, Lambda);
    HTTree res = b.tree()-Ax.tree();
    res.orthogonalize();
    residual   = res.L2normorthogonal()/nrmb;
    #ifdef VERBOSE
        std::cout << "rank1als_sym: Sweep " << 0 << ", leaf " << 0
                  << ", residual = " << residual << "\n";
    #endif
    if (residual<=tol) return 0;

    /* ALS sweeps */
    for (unsigned sweep=1; sweep<=max_sweep; ++sweep) {
        for (unsigned j=1; j<=A.dim(); ++j) {
            /* Compute projection */
            Matrix Pj     = projection(Ax.tree(), x.tree(), j);
            Matrix B      = contract(b.tree(), x.tree(), j);

            /* CG solve on leaf*/
            htucker::DimensionIndex idx(1);
            idx[0]           = j;
            Matrix   xk      = extract(x.tree(), idx);
            T res_cg;
            unsigned it      = cg(A, P, xk, B, Pj, x, j,
                                              Lambda[j-1],
                                              res_cg, tol_cg, maxit_cg);
            #ifdef VERBOSE
                std::cout << "rank1als_sym: CG required " << it
                          << " iterations to reach "
                          << res_cg << "\n";
            #endif

            /* Update x */
            insert(x.tree(), xk, idx);
            if (orthog) x.tree().orthogonalize();

            /* Full residual */
            Ax         = eval(A, x, Lambda, Lambda);
            res        = b.tree()-Ax.tree();
            res.orthogonalize();
            residual   = res.L2normorthogonal()/nrmb;
            #ifdef VERBOSE
                std::cout << "rank1als_sym: Sweep " << sweep << ", leaf " << j
                          << ", residual = " << residual << "\n";
            #endif
            if (residual<=tol) return sweep;
        }
    }

    std::cerr << "rank1als_sym: Reached max sweeps " << max_sweep << "\n";
    return max_sweep;
}


template <typename Optype, typename Prec, typename T, typename Basis>
unsigned
precrank1als_sym(       Sepop<Optype>&                      A,
                        Prec&                               P,
                        HTCoefficients<T, Basis>&           x,
                  const HTCoefficients<T, Basis>&           b,
                  const std::vector<IndexSet<Index1D> >&    Lambda,
                  const flens::GeMatrix
                        <flens::FullStorage
                        <T, cxxblas::ColMajor> >&           scales,
                        T&                                  residual,
                  const bool                                check_res,
                  const bool                                orthog,
                  const bool                                sw,
                  const T                                   balance,
                  const T                                   tol,
                  const unsigned                            max_sweep,
                  const T                                   tol_cg,
                  const unsigned                            maxit_cg,
                  const bool                                verbose_cg)
{
    assert(A.dim()==Lambda.size());
    assert(A.dim()==(unsigned)x.dim());
    assert(x.dim()==b.dim());

    using flens::_;
    typedef typename flens::GeMatrix
                    <flens::FullStorage<T, cxxblas::ColMajor> >     Matrix;
    typedef typename htucker::HTuckerTree<T>                        HTTree;
    typedef HTCoefficients<T, Basis>                                HTCoeff;

    /* Initial residual */
    HTCoeff Ax, xold;
    T nrmb;

    Ax = eval(A, x, Lambda, Lambda);
    if (check_res) {
        HTTree res = b.tree()-Ax.tree();
        res.orthogonalize();
        nrmb       = nrm2(const_cast<HTCoeff&>(b));
        residual   = res.L2normorthogonal()/nrmb;
            std::cout << "precrank1als_sym: Sweep " << 0 << ", leaf " << 0
                      << ", residual = " << residual << "\n";
    }

    /* ALS sweeps */
    int prec = 2;
    T eps;

    Matrix Pj, B;
    bool force_switch = false;
    for (unsigned sweep=1; sweep<=max_sweep; ++sweep) {
        xold = x;

        for (unsigned j=1; j<=A.dim(); ++j) {
            /* Compute projection */
            Pj  = projection(Ax.tree(), x.tree(), j);
            eps = Pj(1, 2)/Pj(1, 1);
            std::cout << "Current  eps is     " << eps << std::endl;

            /* Decide on preconditioning */
            T bm_nonprec = (1.+eps*scales(j, 1)*balance)/(1.+eps*scales(j, 2));
            T bm_prec    = balance*scales(j, 1)/scales(j, 2)*
                           (1.+eps*scales(j, 2)*balance)/
                           (1.+eps*scales(j, 1)*balance);
            if (!force_switch) {
                if (sw && bm_nonprec<bm_prec) prec = 1;
                else                          prec = 2;
            }

            /* Compute rhs */
            B = contract(b.tree(), x.tree(), j);

            /* CG solve on leaf*/
            htucker::DimensionIndex idx(1);
            idx[0]           = j;
            Matrix   xk      = extract(x.tree(), idx);
            T res_cg;
            unsigned it;

            switch (prec) {
                case 1: {
                    #ifdef VERBOSE
                        std::cout << "precrank1als_sym: Solving the non-prec case\n";
                    #endif
                    it = cg(A, xk, B, Pj, x, j,
                            Lambda[j-1], res_cg, tol_cg,
                                                 maxit_cg,
                                                 verbose_cg);
                    if (it==maxit_cg) {
                        prec=2;
                        force_switch = !force_switch;
                    }
                    break;
                }
                case 2: {
                    #ifdef VERBOSE
                        std::cout << "precrank1als_sym: Solving the prec case\n";
                    #endif
                    it  = cg_rank1prec(A, P[j-1], xk, B, Pj, x, j,
                                       Lambda[j-1], res_cg, tol_cg,
                                                            maxit_cg,
                                                            verbose_cg);
                    if (it==maxit_cg) {
                        prec=1;
                        force_switch = !force_switch;
                    }
                    break;
                }
            }

            #ifdef VERBOSE
                std::cout << "precrank1als_sym: Sweep " << sweep << ", leaf "
                          << j << ", CG required " << it
                          << " iterations to reach "
                          << res_cg << "\n";
            #endif

            /* Update x */
            insert(x.tree(), xk, idx);
            if (orthog) x.tree().orthogonalize();

            /* Full residual */
            Ax = eval(A, x, Lambda, Lambda);
            if (check_res) {
                HTTree res = b.tree()-Ax.tree();
                res.orthogonalize();
                residual   = res.L2normorthogonal()/nrmb;
                    std::cout << "precrank1als_sym: Sweep " << sweep << ", leaf " << j
                              << ", residual = " << residual << "\n";
            }
        }

        /* Check for stagnation */
        auto nrmx   = nrm2(x, orthog);
        xold.tree() = x.tree()-xold.tree();
        T diff      = nrm2(xold)/nrmx;
        if (diff<=tol) {
            return sweep;
        }
    }

    std::cerr << "rank1als_sym: Reached max sweeps " << max_sweep << "\n";
    return max_sweep;
}


template <typename Optype, typename T, typename Basis>
unsigned
rank1als_laplace(      Sepop<Optype>&                      A,
                 const std::vector<flens::SyMatrix
                       <flens::FullStorage
                       <T, cxxblas::ColMajor> > >&         Astiff,
                       HTCoefficients<T, Basis>&           x,
                 const HTCoefficients<T, Basis>&           b,
                 const std::vector<IndexSet<Index1D> >&    Lambda,
                 const bool                                check_res,
                 const bool                                orthog,
                 const T                                   tol,
                 const unsigned                            max_sweep,
                 const bool                                verbose)
{
    assert(A.dim()==(unsigned) x.dim());
    assert(x.dim()==b.dim());

    using flens::_;
    typedef typename flens::GeMatrix
                    <flens::FullStorage<T, cxxblas::ColMajor> >     Matrix;
    typedef typename flens::SyMatrix
                    <flens::FullStorage<T, cxxblas::ColMajor> >     SyMatrix;
    typedef typename htucker::HTuckerTree<T>                        HTTree;
    typedef HTCoefficients<T, Basis>                                HTCoeff;

    /* Initial residual */
    HTCoeff Ax, xold;
    T nrmb, residual = -1.;

    Ax = eval(A, x, Lambda, Lambda);
    if (check_res) {
        HTTree res = b.tree()-Ax.tree();
        res.orthogonalize();
        nrmb       = nrm2(const_cast<HTCoeff&>(b));
        residual   = res.L2normorthogonal()/nrmb;
            std::cout << "rank1als_laplace: Sweep " << 0 << ", leaf " << 0
                      << ", residual = " << residual << "\n";
    }

    /* ALS sweeps */
    T eps;
    for (unsigned sweep=1; sweep<=max_sweep; ++sweep) {
        xold = x;

        for (unsigned j=1; j<=A.dim(); ++j) {
            Matrix Pj, B;

            /* Compute projection */
            Pj  = projection(Ax.tree(), x.tree(), j);
            eps = Pj(1, 2)/Pj(1, 1);

            /* Compute rhs */
            B = contract(b.tree(), x.tree(), j);

            /* Solve on leaf */
            htucker::DimensionIndex idx(1);
            idx[0]           = j;

            /* Compute X */
            auto Ap = Astiff[j-1];

            int info;
            Matrix xk;
            if (B.numCols()==1) {
                flens::blas::scal(Pj(1, 2), Ap);
                for (int l=1; l<=Ap.numRows(); ++l) {
                    Ap(l, l) += Pj(1, 1);
                }

                xk   = B;
                info = flens::lapack::posv(Ap, xk);
            } else {
                auto n = B.numRows();
                auto r = B.numCols();

                /* Build Kronecker product (inefficient, replace later) */
                SyMatrix At(n*r, flens::Lower);
                for (int k=1; k<=n*r; ++k) {
                    for (int m=k; m<=n*r; ++m) {
                        unsigned r1 = (m-1)/n + 1;
                        unsigned c1 = (k-1)/n + 1;
                        unsigned r2 = (m-1)%n + 1;
                        unsigned c2 = (k-1)%n + 1;
                        if (r2<c2) {
                            unsigned save = r2;
                            r2 = c2;
                            c2 = save;
                        }

                        At(m, k) = Pj(r1, r+c1)*Ap(r2, c2);
                        if (r2==c2)
                            At(m, k) += Pj(r1, c1);
                    }
                }

                /* Vectorize rhs */
                Matrix vecB(n*r, 1);
                for (int k=1; k<=r; ++k) {
                    vecB(_(n*(k-1)+1, n*k), 1) = B(_, k);
                }

                /* Solve */
                info         = flens::lapack::posv(At, vecB);

                /* Matricisize xk */
                xk.resize(n, r);
                for (int k=1; k<=r; ++k){
                    xk(_, k) = vecB(_(n*(k-1)+1, n*k), 1);
                }
            }
            if (verbose) {
                std::cout << "rank1als_laplace: Sweep " << sweep << ", leaf "
                          << j << ", epsilon = " << eps
                          << ", solver info = " << info << "\n";
            }

            /* Update x */
            insert(x.tree(), xk, idx);
            if (orthog) x.tree().orthogonalize();

            /* Full residual */
            Ax = eval(A, x, Lambda, Lambda);
            if (check_res) {
                HTTree res = b.tree()-Ax.tree();
                res.orthogonalize();
                residual   = res.L2normorthogonal()/nrmb;
                    std::cout << "rank1als_laplace: Sweep " << sweep << ", leaf " << j
                              << ", residual = " << residual << "\n";
            }
        }

        /* Check for stagnation */
        auto nrmx   = nrm2(x, orthog);
        xold.tree() = x.tree()-xold.tree();
        T diff      = nrm2(xold)/nrmx;
        if (diff<=tol) {
            return sweep;
        }
    }

    std::cerr << "rank1als_laplace: Reached max sweeps " << max_sweep << "\n";
    return max_sweep;
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_RANK1ALS_TCC
