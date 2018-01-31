#include <iostream>
#include <lawa/lawa.h>
#include <vector>
#include <cmath>
#include <cstdlib>

typedef double                                  T;
typedef flens::DenseVector<flens::Array<T> >    DenseVector;
typedef flens::DenseVector<flens::Array<int> >  IDV;
typedef flens::GeMatrix<flens::FullStorage
        <T, cxxblas::ColMajor> >                GeMat;
typedef flens::GeMatrix<flens::FullStorage
        <int, cxxblas::ColMajor> >              MatInt;
typedef lawa::Function<T>                       Function;
typedef lawa::Basis<T, lawa::Orthogonal,
        lawa::Interval, lawa::Multi>            Basis;
typedef lawa::IndexD                            IndexD;
typedef lawa::Index1D                           Index1D;
typedef lawa::Index2D                           Index2D;
typedef lawa::IndexSet<Index1D>                 IndexSet;
typedef std::vector<IndexSet>                   IndexSetVec;
typedef lawa::RHSWithPeaks1D<T, Basis>          Integral1D;
typedef lawa::Coefficients<
        lawa::Lexicographical, T,
        lawa::Index1D>                          Coeff1D;
typedef lawa::Coefficients<
        lawa::Lexicographical, T,
        lawa::Index2D>                          Coeff2D;
typedef lawa::SepCoefficients<
        lawa::Lexicographical, T, Index1D>      SepCoeff;


typedef     lawa::AdaptiveLaplaceOperator1D<T,
            lawa::Orthogonal,
            lawa::Interval,
            lawa::Multi>             Laplace1D;
typedef     lawa::AdaptiveLaplaceOperator1D<T,
            lawa::Orthogonal,
            lawa::Interval,
            lawa::MultiRefinement>   RefLaplace1D;
typedef     lawa::AdaptiveIdentityOperator1D<T,
            lawa::Orthogonal,
            lawa::Interval,
            lawa::Multi>             Identity1D;
typedef     lawa::AdaptiveIdentityOperator1D<T,
            lawa::Orthogonal,
            lawa::Interval,
            lawa::MultiRefinement>   RefIdentity1D;

typedef     lawa::LocalOperator1D<Basis,
            Basis, RefLaplace1D, Laplace1D>     LOp_Lapl1D;
typedef     lawa::LocalOperator1D<Basis,
            Basis, RefIdentity1D,
            Identity1D>                         LOp_Id1D;

typedef     lawa::Sepop<LOp_Lapl1D>             Sepop;
double
one(double)
{
    return 1.;
}


double
zero(double)
{
    return 0.;
}


unsigned
richardson(Sepop&                          A,
           lawa::Sepdiagscal<Basis>&       S,
           lawa::HTCoefficients<T, Basis>& x,
           lawa::HTCoefficients<T, Basis>& b,
           IndexSetVec&                    Lambda,
           T                               omega,
           T                               tol,
           unsigned                        maxit)
{
    auto r     = b;
    T residual = nrm2(r);
    T delta    = 1e-01;

    std::cout << "richardson: Iteration 0, residual "
              << residual << std::endl;
    if (residual<=tol) return 0;

    for (unsigned k=1; k<=maxit; ++k) {
        scal(omega, r);
        T trunc  = delta*omega*residual;
        if (k==1) x   = r;
        else x.tree() = add_truncate(x.tree(), r.tree(), trunc);

        r        = evaleff2(A, S, x, Lambda, Lambda, trunc/2.);
        scal(-1., r);
        r.tree() = add_truncate(b.tree(), r.tree(), trunc/2.);
        residual = nrm2(r);

        std::cout << "richardson: Iteration " << k
                  << ", residual " << residual << std::endl;
    }

    return maxit;
}


unsigned
htawgm3(lawa::Sepdiagscal<Basis>&       S,
        Sepop&                          A,
        lawa::HTCoefficients<T, Basis>& u0,
        lawa::SeparableRHSD<T, Basis>   f,
        IndexSetVec&                    Lambda0,
        T                               omega0,
        T                               omega1,
        T                               omega2,
        T                               omega3,
        T                               omega4,
        T                               omega5,
        T                               lambda_max,
        T                               lambda_min,
        T                               alpha,
        T                               tol,
        unsigned                        I,
        unsigned                        M,
        unsigned                        K,
        T&                              residual,
        bool                            verbose)
{
    // Input check
    assert(S.dim()==A.dim());
    assert(S.dim()==(unsigned) u0.dim());
    assert(S.dim()==f.dim());
    assert(S.dim()==Lambda0.size());

    // Initial residual and rhs
    residual = omega0;
    T fdelta = omega0;
    if (verbose) {
        std::cout << "htawgm: Iteration 0, residual " << residual
                  << std::endl;
    }
    if ((1.+omega1)*residual<=tol) return 0;

    // Timing tools
    auto start   = std::chrono::system_clock::now();
    auto end     = start;
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start);

    // Outer iterations
    bool zero = true;
    IndexSetVec sweep = Lambda0;
    IndexSetVec total = Lambda0;
    lawa::HTCoefficients<T, Basis> rtree(u0.dim(),
                                         u0.basis(),
                                         u0.map());
    lawa::HTCoefficients<T, Basis> f_Lambda(u0.dim(),
                                            u0.basis(),
                                            u0.map());
    SepCoeff fcp(f.rank(), f.dim());
    for (unsigned k=0; k<=K; ++k) {
        // Inner iterations
        for (unsigned m=1; m<=M; ++m) {
            // Index set size and rank
            if (verbose) {
                std::cout << "htawgm: Iteration (" << k
                          << ", " << m << ")\n";
                for (unsigned j=0; j<Lambda0.size(); ++j) {
                    unsigned jmax = 0;
                    for (const auto& it : Lambda0[j]) {
                        unsigned l = it.j;
                        if (it.xtype==lawa::XWavelet) ++l;
                        jmax = std::max(jmax, l);
                    }

                    std::cout << "        d = " << j+1
                              << ", size = " << Lambda0[j].size()
                              << ", max level = " << jmax << std::endl;
                }

                std::cout << "        rank = " << u0.tree().max_rank()
                          << std::endl;
            }
            // Fix rhs for pcg
            start   = std::chrono::system_clock::now();

            genCoefficients(fcp, f, Lambda0);
            set(f_Lambda, fcp);

            // Determine scaling precision
            T eps_km = omega2*residual;
            T eta_km = (1.-S.eps())*eps_km/
                       (3.*(2.*fdelta+2.*lambda_max*nrm2(u0)));
            S.set_nu(eta_km);
            T iscale = compIndexscale(u0.basis(), Lambda0, S.order());
            S.set_iscale(iscale);
            S.comp_n();
            f_Lambda = applyScale(S, f_Lambda, Lambda0, 0.5*eps_km);

            end     = std::chrono::system_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::seconds>
                      (end-start);
            if (verbose) {
                    std::cout << "htawgm: Iteration (" << k
                          << ", " << m << ")\n";
                    std::cout << "        Set up rhs took "
                              << elapsed.count() << " secs\n";
            }

            if (verbose) {
                std::cout << "htawgm: Iteration (" << k
                          << ", " << m << ")\n";
                std::cout << "        eps_km  = " << eps_km << "\n";
                std::cout << "        eta_km  = " << eta_km << "\n";
                std::cout << "        rank(S) = "
                          << S.n()+S.nplus()+1 << "\n";
            }

            // PCG iterations
            start   = std::chrono::system_clock::now();

            if (verbose) {
                 std::cout << "htawgm: Iteration (" << k
                           << ", " << m << ")\n";
                 std::cout << "\n*** PCG begin ***\n\n";
            }
            T        res_pcg;
            T        delta  = 1e-01;
            unsigned pcgit = galerkin_pcg2(A, S, u0, f_Lambda, Lambda0,
                                          res_pcg,
                                          zero, eps_km, I, delta, 0.5*eps_km);
            zero            = false;
            if (verbose) {
                 std::cout << "htawgm: pcg required " << pcgit
                           << " iterations for r<="   << res_pcg
                           << std::endl;
                 std::cout << "\n*** PCG end ***\n\n";
            }

            end     = std::chrono::system_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::seconds>
                      (end-start);
            if (verbose) {
                    std::cout << "htawgm: Iteration (" << k
                          << ", " << m << ")\n";
                    std::cout << "        PCG took "
                              << elapsed.count() << " secs\n";
            }
            // Evaluate residual
            start  = std::chrono::system_clock::now();

            eps_km = omega1*eps_km;
            eta_km = (1.-S.eps())*eps_km/
                     (3.*(2.*fdelta+2.*lambda_max*nrm2(u0)));
            S.set_nu(eta_km);

            sweep    = presidual2(A, S, u0, f_Lambda, fcp, rtree, f,
                                  Lambda0, sweep, total,
                                  eps_km);
            residual = nrm2(rtree);

            end      = std::chrono::system_clock::now();
            elapsed  = std::chrono::duration_cast<std::chrono::seconds>
                      (end-start);

            if (verbose) {
                    std::cout << "htawgm: Iteration (" << k
                          << ", " << m << ")\n";
                    std::cout << "        RES took "
                              << elapsed.count() << " secs\n";
            }
            if (verbose) {
                 std::cout << "htawgm: Iteration (" << k
                           << ", " << m << "), residual = "
                           << residual << "\n";
            }
            if ((1.+omega1)*residual<=tol) return k;

            // Check stopping criterion
            if ((1.+omega1)*residual<=omega3*omega0) break;

            // Extend index set
            start   = std::chrono::system_clock::now();

            sweep   = bulk(alpha, residual, rtree, Lambda0, sweep);
            extend(u0, Lambda0);

            end     = std::chrono::system_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::seconds>
                      (end-start);

            if (verbose) {
                    std::cout << "htawgm: Iteration (" << k
                          << ", " << m << ")\n";
                    std::cout << "        EXPAND took "
                              << elapsed.count() << " secs\n";
            }
        }

        // Truncate and coarsen
        T epst = omega4*omega0/lambda_min;
        T epsc = omega5*omega0/lambda_min;

        if (verbose) {
            std::cout << "htawgm: Iteration "    << k        << "\n";
            std::cout << "        eps trunc  = " << epst     << "\n";
            std::cout << "        eps coarse = " << epsc     << "\n";
            std::cout << "        nrm(u)     = " << nrm2(u0) << "\n";
        }

        start   = std::chrono::system_clock::now();

        u0.truncate(epst);
        Lambda0 = coarsen(u0, Lambda0, epsc);
        sweep   = Lambda0;
        total   = Lambda0;

        end     = std::chrono::system_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::seconds>
                  (end-start);

        if (verbose) {
                std::cout << "htawgm: Iteration " << k << "\n";
                std::cout << "        C(T(...)) took "
                          << elapsed.count() << " secs\n";
        }

        (void) presidual2(A, S, u0, f_Lambda, fcp, rtree, f,
                          Lambda0, sweep, total,
                          residual);
        residual = nrm2(rtree);
        if (verbose) {
             std::cout << "htawgm: Iteration " << k+1
                       << ", residual = "
                       << residual << "\n";
        }

        // Update omega0
        omega0 *= (omega3+omega4+omega5);
    }

    std::cerr << "htawgm: Warning! Max iterations reached!\n";
    return K;
}

T
alpha_(T x)
{
    return std::pow(std::log(1.+std::exp(x)), 2.);
}

T
omega_(T x)
{
    return 2./std::sqrt(M_PI)*(1./(1.+std::exp(-1.*x)));
}

T
appScale(T t, T delta, T eta, T scale, int dim = 1)
{
    T h      = std::pow(M_PI, 2.)/(5.*(std::abs(std::log(delta/2.))+4.));
    T max    = std::max(4./std::sqrt(M_PI), std::sqrt(std::abs
               (std::log(delta/2.))));
    T min    = std::min(delta/2., eta);
    int np   = std::ceil(max/h);
    int n    = std::ceil((1./h)*(std::log(2.*std::pow(M_PI, -0.5))+
                         +std::abs(std::log(min))
                         +0.5*std::log(scale)));


    T varphi = 0.;
    for (int i=-n; i<=np; ++i) {
        varphi += std::pow(omega_(h*(T) i)*h, 1./(T) dim)
                 *std::exp(-alpha_(h*(T) i)*t);
    }

    return varphi;
}


T
factor(T t, T omega_min, T delta, int k, int dim = 1)
{
    T h      = std::pow(M_PI, 2.)/(5.*(std::abs(std::log(delta/2.))+4.));
    T varphi = std::pow(omega_(h*(T) k)*h/omega_min, 1./(T) dim)
                 *std::exp(-alpha_(h*(T) k)*t);

    return varphi;
}


int
main(int argc, char* argv[])
{
    if (argc!=4) {
        std::cerr << "Usage: " << argv[0] << " dim level rich/awgm(0/1)\n";
        return 1;
    }

    int dim  = atoi(argv[1]);
    int lev  = atoi(argv[2]);
    int awgm = atoi(argv[3]);
    if (dim<=0) {
        std::cerr << "Dimension must be a positive integer\n";
        return 1;
    }
    if (!(awgm==1 || awgm==0)) {
        std::cerr << "Third argument must be 0 or 1\n";
        return 1;
    }

    int rank = 1;

    DenseVector     sings;
    GeMat           deltas(3,2);

    std::vector<GeMat> _deltas;
    Function        onef(one, sings);
    Basis           basis(4);
    basis.template enforceBoundaryCondition<lawa::DirichletBC>();

    if (lev<basis.j0) {
        std::cerr << "Level too small\n";
        return 1;
    }

    IndexSet        indexset;
    std::vector<Function> fvec;

    for (int i=1; i<=dim; ++i) {
        fvec.push_back(onef);
    }

    SepCoeff        coeffs(rank, dim);
    IndexSetVec     indexsetvec(dim);
    lawa::SeparableFunctionD<T> F(fvec, rank, dim);
    MatInt                      derivs(rank, dim);
    for (int i=1; i<=rank; ++i) {
        for (int j=1; j<=dim; ++j) {
            derivs(i,j) = 0;
            _deltas.push_back(deltas);
        }
    }

    lawa::SeparableRHSD<T, Basis>   Fint(basis, F, _deltas, derivs);

    getFullIndexSet(basis, indexset,  lev);

    std::cout << "The initial index set size is " << indexset.size()
              << "\n\n";

    for (int l=0; (unsigned)l<indexsetvec.size(); ++l) {
        indexsetvec[l]   = indexset;
    }

    /* Map */
    lawa::Mapwavind<Index1D> map(dim);

    lawa::HTCoefficients<T, Basis>    f(dim, basis, map);
    lawa::HTCoefficients<T, Basis>    u(dim, basis, map);
    lawa::HTCoefficients<T, Basis>    r(dim, basis, map);

    Laplace1D       LaplaceBil(basis);
    RefLaplace1D    RefLaplaceBil(basis.refinementbasis);
    Identity1D      IdentityBil(basis);
    RefIdentity1D   RefIdentityBil(basis.refinementbasis);
    LOp_Lapl1D      lapl(basis, basis, RefLaplaceBil, LaplaceBil);

    rndinit(u, indexsetvec, 2, 0.5);
    Sepop A(lapl, dim, dim);

    lawa::Sepdiagscal<Basis>    S(dim, basis);

    lawa::HTRICH_Params  params;
    params.omega      = 1e-03;
    params.cA         = 1.;
    params.eps0       = 5e-02;
    params.rho        = 0.92;
    params.maxit_inner= 250;
    params.beta1      = 1e-01;
    params.beta2      = 1e-01;
    params.maxit_rich = 100;

    double alpha      = 1e-02;
    double kappaP     = std::sqrt(2.*dim-3.);
    double kappaC     = std::sqrt(dim);
    params.kappa1     = 1./(1.+(1.+alpha)*(kappaP+kappaC+kappaP*kappaC));
    params.kappa2     = (1.+alpha)*kappaP*params.kappa1;
    params.kappa3     = kappaC*(kappaP+1.)*(1.+alpha)*params.kappa1;

    lawa::HTAWGM_Params params2;
    params2.nrmA       = 100.;
    params2.maxit_pcg  = 50;
    params2.maxit_awgm = 300;
    params2.tol_awgm   = 1e-08;
    params2.delta1_pcg = .1;
    params2.alpha      = 0.5;
    params2.gamma0     = 1e-01;
    params2.gamma1     = 8e-02;
    params2.gammait    = 9;

//    std::cout << "HTRICH params =\n";
//    std::cout << params << std::endl;
//
//    std::cout << "HTAWGM params =\n";
//    std::cout << params2 << std::endl;

    unsigned its = 0;
    double   res = 0;

    genCoefficients(coeffs, Fint, indexsetvec);
    set(f, coeffs);
    //rndinit(f, indexsetvec, 5, 1.);
    T delta  = 1e-01;
    T eta    = 1e-01;
    setScaling(S, delta);
    S.set_nu(eta);

    /*
    if (!awgm) {
        its = htrich(A, S, u, Fint, indexsetvec, res, params);
    } else {
        its = htawgm2(A, S, u, Fint, indexsetvec, res, params2);
    }*/

    T lambda_max = 2.;
    T lambda_min = 0.1;
    T kA         = lambda_max/lambda_min;
    T comprho    = 1.01;
    T outrho     = 0.1;
    alpha        = 0.8;

    T omega0     = 0.65;
    T omega1     = 0.25;
    T omega2     = 0.1;
    T omega3     = outrho/(1.+comprho*std::sqrt(2.*dim-3)+
                              comprho*std::sqrt((T) dim)*
                              (1.+std::sqrt(2.*dim-3)));
    T omega4     = comprho*std::sqrt(2.*dim-3)*omega3;
    T omega5     = comprho*std::sqrt((T) dim)*(1.+std::sqrt(2.*dim-3))*omega3;
    T tol        = 1e-08;

    T inrho      = std::sqrt(1.-std::pow((alpha-omega1)/(1.+omega1), 2.)/kA+
                             std::pow(omega2/(1.-omega1), 2.)*kA);

    unsigned I = 50;
    unsigned M = 50;
    unsigned K = 50;

    bool verbose = true;

    std::cout << "*** HTAWGM Parameters ***\n";
    std::cout << "lambda_max = " << lambda_max << std::endl;
    std::cout << "lambda_min = " << lambda_min << std::endl;
    std::cout << "kA         = " << kA << std::endl;
    std::cout << "alpha      = " << alpha << std::endl;
    std::cout << "omega0     = " << omega0 << std::endl;
    std::cout << "omega1     = " << omega1 << std::endl;
    std::cout << "omega2     = " << omega2 << std::endl;
    std::cout << "omega3     = " << omega3 << std::endl;
    std::cout << "omega4     = " << omega4 << std::endl;
    std::cout << "omega5     = " << omega5 << std::endl;
    std::cout << "vartheta   = " << inrho  << std::endl;
    std::cout << "*** ----------------- ***\n";

 // Compare to exact scaling
 //   std::cout << "rank f = " << f.tree().max_rank() << std::endl;
 //   std::cout << "nrm(f) = " << nrm2(f) << std::endl;
 //   f.tree().print_w_UorB();
 //   auto ist   = applyScale(S, f, indexsetvec, 1e-08/nrm2(f));
 //   S.set_nu(1e-07);
 //   auto soll2 = applyScale(S, f, indexsetvec, 1e-08/nrm2(f));
 //   lawa::Coefficients<lawa::Lexicographical, T, lawa::Index2D> soll;
 //   T omega_min = compOmegamin2(S.basis(), S.dim(), S.order());

 //   for (auto& mu1 : indexsetvec[0]) {
 //       for (auto& mu2 : indexsetvec[1]) {
 //           lawa::IndexD  indexd(2);
 //           lawa::Index2D index(mu1, mu2);
 //           indexd(1)   = mu1;
 //           indexd(2)   = mu2;
 //           int l1      = mu1.j;
 //           int l2      = mu2.j;
 //           if (mu1.xtype==lawa::XWavelet) ++l1;
 //           if (mu2.xtype==lawa::XWavelet) ++l2;

 //           T scale1 = std::pow(2., 2.*(T) l1);
 //           T scale2 = std::pow(2., 2.*(T) l2);
 //           T scale  = scale1+scale2;

 //           soll[index] = f(indexd)/std::sqrt(scale);
 //       }
 //   }

 //   T err  = 0.;
 //   T err2 = 0.;
 //   soll2.tree() = ist.tree()-soll2.tree();
 //   soll2.orthogonalize();
 //   for (auto& mu : soll) {
 //       lawa::IndexD  indexd(2);
 //       indexd(1) = mu.first.index1;
 //       indexd(2) = mu.first.index2;
 //       int l1    = mu.first.index1.j;
 //       int l2    = mu.first.index2.j;
 //       if (mu.first.index1.xtype==lawa::XWavelet) ++l1;
 //       if (mu.first.index1.xtype==lawa::XWavelet) ++l2;

 //       T scale     = std::pow(2., 2.*(T) l1);
 //         scale    += std::pow(2., 2.*(T) l2);

 //       err  += scale*std::pow(mu.second-ist(indexd), 2.);
 //       err2 += scale*std::pow(soll2(indexd), 2.);
 //   }
 //   std::cout << "*\n*\n*\nErr1 = "
 //             << std::sqrt(err)/nrm2(f) << std::endl;
 //   std::cout << "Err2 = "
 //             << std::sqrt(err2)/nrm2(f) << std::endl;


 //   T t      = (std::pow(2., 2.)+std::pow(2., 2.*2.));
 //   T c      = std::sqrt(omega_min);

 //   T varphi = appScale(t/(c*c), delta, eta, t/(c*c));
 //   err = std::sqrt(t)*std::abs(1./(std::sqrt(t))-varphi/c);
 //   std::cout << "Err3 = " << err << std::endl;

 //   exit(1);

    its = htawgm3(S, A, u, Fint, indexsetvec,
                  omega0,
                  omega1,
                  omega2,
                  omega3,
                  omega4,
                  omega5,
                  lambda_max,
                  lambda_min,
                  alpha,
                  tol,
                  I,
                  M,
                  K,
                  res,
                  verbose);

    std::cout << "Solver took " << its << " iterations to reach "
              << res << " accuracy" << std::endl;

    return 0;
}
