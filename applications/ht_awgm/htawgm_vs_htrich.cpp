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
typedef lawa::Coefficients<
        lawa::Lexicographical, T,
        lawa::IndexD>                           CoeffD;
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

lawa::NoPreconditioner<T, Index1D>          pnull;

void
saveCoeffVector1D(const Coeff1D &coeff, const Basis &basis, const char* filename)
{
    typedef typename Coeff1D::const_iterator const_coeff_it;

    std::ofstream data(filename);

    if(!data.good()){
        std::cerr << "File "
                  << filename
                  << " could not be opened for writing" << std::endl;
        exit(1);
    }
    data.precision(40);
    data << "# Center_x Value Xtype j k" << std::endl;

    for (const_coeff_it it = coeff.begin(); it != coeff.end(); ++it) {
      FLENS_DEFAULT_INDEXTYPE j=(*it).first.j,
                              k=(*it).first.k;
      lawa::XType type=(*it).first.xtype;

      //center of the support
      double x = 0.5*(basis.generator(type).support(j,k).l2 +
                      basis.generator(type).support(j,k).l1);

      data << x << " " << std::scientific <<  (*it).second << " "
           << type << " " << j << " " << k << std::endl;
    }
    data.close();
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
    lawa::HTCoefficients<T, Basis> f_Lambda(u0.dim(),
                                            u0.basis(),
                                            u0.map());
    SepCoeff fcp(f.rank(), f.dim());

    genCoefficients(fcp, f, Lambda0);
    set(f_Lambda, fcp);
    S.assemble(Lambda0);
    f_Lambda     = applyScale(S, f_Lambda, Lambda0, tol);

    T fdelta     = nrm2(f_Lambda);
    residual     = fdelta;
    T omega0     = residual;
    auto Ttotal0 = std::chrono::system_clock::now();
    auto Ttotal1 = std::chrono::system_clock::now();
    auto dt      = std::chrono::duration_cast<std::chrono::seconds>
                   (Ttotal1-Ttotal0);
    if (verbose) {
        std::cout << "htawgm: Iteration 0, residual " << residual
                  << std::endl;
        std::cout << "        Current time: " << dt.count()
                  << " secs\n";
    }
    if ((1.+omega1)*residual<=tol) return 0;

    // Timing tools
    auto start   = std::chrono::system_clock::now();
    auto end     = start;
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start);

    // Outer iterations
    bool zero = true;
    IndexSetVec Lambda = Lambda0;
    IndexSetVec sweep  = Lambda0;
    IndexSetVec total  = Lambda0;
    lawa::HTCoefficients<T, Basis> rtree(u0.dim(),
                                         u0.basis(),
                                         u0.map());
    T err0;
    for (unsigned k=0; k<=K; ++k) {
        // Inner iterations
        for (unsigned m=1; m<=M; ++m) {
            // Index set size and rank
            if (verbose) {
                std::cout << "htawgm: Iteration (" << k
                          << ", " << m << ")\n";
                unsigned sum = 0;
                for (unsigned j=0; j<Lambda.size(); ++j) {
                    unsigned jmax = 0;
                    for (const auto& it : Lambda[j]) {
                        unsigned l = it.j;
                        if (it.xtype==lawa::XWavelet) ++l;
                        jmax = std::max(jmax, l);
                    }

                    sum += Lambda[j].size();
                    std::cout << "        d = " << j+1
                              << ", size = " << Lambda[j].size()
                              << ", max level = " << jmax << std::endl;
                }

                std::cout << "        total size = " << sum << std::endl;
            }
            // Fix rhs for pcg
            start   = std::chrono::system_clock::now();

            genCoefficients(fcp, f, Lambda);
            set(f_Lambda, fcp);

            // Determine scaling precision
            T lower  = 10.;
            T eps_km = omega2*residual;
            T eta_km = (1.-S.eps())*residual*lower/
                       (3.*(2.*fdelta+2.*lambda_max*nrm2(u0)));
            //S.set_nu(eta_km);
            S.assemble(Lambda);
            f_Lambda = applyScale(S, f_Lambda, Lambda, eps_km);

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
            unsigned pcgit  = galerkin_pcg2(A, S, u0, f_Lambda, Lambda,
                                            res_pcg,
                                            zero, eps_km, I, delta, residual);
            zero            = false;
            if (k==0 && m==1) {
                err0 = nrm2(u0);
                std::cout << "err0 = " << err0 << std::endl;
            }

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
            start    = std::chrono::system_clock::now();

            sweep    = presidual2(A, S, u0, f_Lambda, fcp, rtree, f,
                                  Lambda, sweep, total,
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
                 Ttotal1 = std::chrono::system_clock::now();
                 dt      = std::chrono::duration_cast<std::chrono::seconds>
                          (Ttotal1-Ttotal0);
                 std::cout << "        Current time: " << dt.count()
                          << " secs\n";

            }
            if (verbose) {
                std::cout << "htawgm: rank = " << u0.tree().max_rank()
                          << std::endl;
            }

            if ((1.+omega1)*residual<=tol) return k;

            // Check stopping criterion
            if ((1.+omega1)*residual<=omega3*err0) break;

            // Extend index set
            start   = std::chrono::system_clock::now();

            sweep   = bulk(alpha, residual, rtree, Lambda, sweep);
            extend(u0, Lambda);

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

        Lambda0 = Lambda;
        return k;

        // Truncate and coarsen
        T epst = omega4*err0;
        T epsc = omega5*err0;

        if (verbose) {
            std::cout << "htawgm: Iteration "    << k        << "\n";
            std::cout << "        eps trunc  = " << epst     << "\n";
            std::cout << "        eps coarse = " << epsc     << "\n";
            std::cout << "        nrm(u)     = " << nrm2(u0) << "\n";
        }

        start   = std::chrono::system_clock::now();

        u0.truncate(epst);
        Lambda  = coarsen(u0, Lambda, epsc);
        Lambda  = unify(Lambda0, Lambda);
        restrict(u0, Lambda);
        sweep   = Lambda;
        total   = Lambda;

        end     = std::chrono::system_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::seconds>
                  (end-start);

        if (verbose) {
                std::cout << "htawgm: Iteration " << k << "\n";
                std::cout << "        C(T(...)) took "
                          << elapsed.count() << " secs\n";
        }
        if (verbose) {
                std::cout << "htawgm: rank = " << u0.tree().max_rank()
                          << std::endl;
        }

        genCoefficients(fcp, f, Lambda);
        omega0 *= omega3+omega4+omega5;

        // Determine scaling precision
        T eps_km = residual;
        T eta_km = (1.-S.eps())*residual*10/
                   (3.*(2.*fdelta+2.*lambda_max*nrm2(u0)));
        //S.set_nu(eta_km);
        S.assemble(Lambda);
        (void) presidual2(A, S, u0, f_Lambda, fcp, rtree, f,
                          Lambda, sweep, total,
                          eps_km);
        residual = nrm2(rtree);

        sweep   = Lambda;
        total   = Lambda;

        if (verbose) {
             std::cout << "htawgm: Iteration " << k+1
                       << ", omega0 = " << omega0
                       << ", residual = "
                       << residual << "\n";
             Ttotal1 = std::chrono::system_clock::now();
             dt      = std::chrono::duration_cast<std::chrono::seconds>
                      (Ttotal1-Ttotal0);
             std::cout << "        Current time: " << dt.count()
                      << " secs\n";

        }

        // Update err0
        err0 *= omega3+omega4+omega5;
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

    T sp = 0.5;
    lawa::HTCoefficients<T, Basis>    f(dim, sp, basis, map);
    lawa::HTCoefficients<T, Basis>    u(dim, sp, basis, map);
    lawa::HTCoefficients<T, Basis>    r(dim, sp, basis, map);

    Laplace1D       LaplaceBil(basis);
    RefLaplace1D    RefLaplaceBil(basis.refinementbasis);
    Identity1D      IdentityBil(basis);
    RefIdentity1D   RefIdentityBil(basis.refinementbasis);
    LOp_Lapl1D      lapl(basis, basis, RefLaplaceBil, LaplaceBil);

    rndinit(u, indexsetvec, 1, 0.);
    Sepop A(lapl, dim, dim);

    lawa::Sepdiagscal<Basis>    S(dim, basis, map);

//    lawa::HTRICH_Params  params;
//    params.omega      = 1e-03;
//    params.cA         = 1.;
//    params.eps0       = 5e-02;
//    params.rho        = 0.92;
//    params.maxit_inner= 250;
//    params.beta1      = 1e-01;
//    params.beta2      = 1e-01;
//    params.maxit_rich = 100;
//
//    double alpha      = 1e-02;
//    double kappaP     = std::sqrt(2.*dim-3.);
//    double kappaC     = std::sqrt(dim);
//    params.kappa1     = 1./(1.+(1.+alpha)*(kappaP+kappaC+kappaP*kappaC));
//    params.kappa2     = (1.+alpha)*kappaP*params.kappa1;
//    params.kappa3     = kappaC*(kappaP+1.)*(1.+alpha)*params.kappa1;

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

    unsigned its = 0;
    double   res = 0;

    genCoefficients(coeffs, Fint, indexsetvec);
    set(f, coeffs);
    T delta  = 0.5;
    T eta    = 0.1;
    setScaling(S, delta);
    S.set_nu(eta);

    /*
    if (!awgm) {
        its = htrich(A, S, u, Fint, indexsetvec, res, params);
    } else {
        its = htawgm2(A, S, u, Fint, indexsetvec, res, params2);
    }*/

    T lambda_max = 3.;
    T lambda_min = 0.3;
    T kA         = lambda_max/lambda_min;
    T comprho    = 1.01;
    T outrho     = 0.1;
    T alpha      = 0.9;

    T omega1     = 0.1;
    T omega2     = 0.5;
    T omega3     = outrho/(1.+comprho*std::sqrt(2.*dim-3)+
                              comprho*std::sqrt((T) dim)*
                              (1.+std::sqrt(2.*dim-3)));
    T omega4     = comprho*std::sqrt(2.*dim-3)*omega3;
    T omega5     = comprho*std::sqrt((T) dim)*(1.+std::sqrt(2.*dim-3))*omega3;
    T tol        = 1e-06;

    unsigned I = 10;
//    unsigned M = std::ceil(std::abs(std::log(omega3/std::sqrt(kA))/
//                           std::log(inrho)));
    unsigned M = 15;
//    unsigned K = std::ceil(std::abs(
//                 std::log(1./(tol*kA*omega3*omega0*(1.+omega1))*(1.-omega1))/
//                 std::log(omega3+omega4+omega5)));
    unsigned K = 0;

    bool verbose = true;

    std::cout << "*** HTAWGM Parameters ***\n";
    std::cout << "lambda_max = " << lambda_max << std::endl;
    std::cout << "lambda_min = " << lambda_min << std::endl;
    std::cout << "kA         = " << kA << std::endl;
    std::cout << "alpha      = " << alpha << std::endl;
    std::cout << "omega1     = " << omega1 << std::endl;
    std::cout << "omega2     = " << omega2 << std::endl;
    std::cout << "omega3     = " << omega3 << std::endl;
    std::cout << "omega4     = " << omega4 << std::endl;
    std::cout << "omega5     = " << omega5 << std::endl;
    std::cout << "I          = " << I      << std::endl;
    std::cout << "M          = " << M      << std::endl;
    std::cout << "K          = " << K      << std::endl;
    std::cout << "*** ----------------- ***\n";

//    T trunc = 1e-02*nrm2(u);
//
//    auto t0 = std::chrono::high_resolution_clock::now();
//    u.truncate(trunc);
//    auto t1 = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<T> dt = t1-t0;
//    std::cout << "trunc took " << dt.count() << " secs\n";
//
//
//    exit(1);
//    its = htawgm3(S, A, u, Fint, indexsetvec,
//                  omega1,
//                  omega2,
//                  omega3,
//                  omega4,
//                  omega5,
//                  lambda_max,
//                  lambda_min,
//                  alpha,
//                  tol,
//                  I,
//                  M,
//                  K,
//                  res,
//                  verbose);

    std::cout << "Solver took " << its << " iterations to reach "
              << res << " accuracy" << std::endl;

    htucker::DimensionIndex idx(1);
    idx[0]    = 1;
    auto outp = extract(u, indexsetvec[0], idx);

    saveCoeffVector1D(outp(1, 1), u.basis(), "support_d1.dat");
    exit(1);
    for (unsigned k=1; k<=outp.rank(); ++k) {
        std::string name = "plot_";
        name            += std::to_string(k);
        name            += ".dat";
        plot(basis, outp(k, 1), pnull, zero, zero, 0., 1.02, 1./100.,
             name.c_str());
    }

    return 0;
}
