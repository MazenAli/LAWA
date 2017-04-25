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
typedef lawa::IndexSet<Index1D>                 IndexSet;
typedef std::vector<IndexSet>                   IndexSetVec;
typedef lawa::RHSWithPeaks1D<T, Basis>          Integral1D;
typedef lawa::Coefficients<
        lawa::Lexicographical, T,
        lawa::Index1D>                          Coeff1D;
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

typedef     lawa::Sepop<lawa::AbstractLocalOperator1D<T, Basis>>
                                                Sepop;
double
one(double)
{
    return 1.;
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
    IndexSet        indexset2;
    Coeff1D         coeff;

    std::vector<Function> fvec;

    for (int i=1; i<=dim; ++i) {
        fvec.push_back(onef);
    }

    SepCoeff        coeffs(rank, dim);
    IndexSetVec     indexsetvec(dim);
    IndexSetVec     indexsetvec2(dim);
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
    getFullIndexSet(basis, indexset2,  lev+2);

    std::cout << "The initial index set size is " << indexset.size()
              << "\n\n";

    for (int l=0; (unsigned)l<indexsetvec.size(); ++l) {
        indexsetvec[l] = indexset;
        indexsetvec2[l] = indexset2;
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

    Sepop A(lapl, dim, dim);

    lawa::Sepdiagscal<Basis>    S(dim, basis);

    lawa::HTRICH_Params  params;
    params.omega      = 1e-03;
    params.cA         = 1.;
    params.eps0       = 1e-02;
    params.rho        = 0.93;
    params.maxit_inner= 60;
    params.beta1      = 1e-02;
    params.beta2      = 1e-02;
    params.maxit_rich = 100;

    double alpha      = 1e-02;
    double kappaP     = std::sqrt(2.*dim-3.);
    double kappaC     = std::sqrt(dim);
    params.kappa1     = 1e-01/(1.+(1.+alpha)*(kappaP+kappaC+kappaP*kappaC));
    params.kappa2     = (1.+alpha)*kappaP*params.kappa1;
    params.kappa3     = 1e-01*kappaC*(kappaP+1.)*(1.+alpha)*params.kappa1;

    lawa::HTAWGM_Params params2;
    params2.maxit_pcg  = 50;
    params2.maxit_awgm = 100;
    params2.tol_awgm   = 1e-08;
    params2.delta1_pcg = 1e-01;
    params2.delta2_pcg = 1e-01;
    params2.delta3_pcg = 1e-01;
    params2.alpha      = 0.95;
    params2.gamma0     = 1e-02;
    params2.gamma1     = 0.25;
    params2.gammait    = 9;

    std::cout << "HTRICH params =\n";
    std::cout << params << std::endl;

    std::cout << "HTAWGM params =\n";
    std::cout << params2 << std::endl;

    unsigned its = 0;
    double   res = 0;

    genCoefficients(coeffs, Fint, indexsetvec);
    set(f, coeffs);
    setScaling(S, 5e-01);
    S.set_nu(1e-01);

    if (!awgm) {
        its = htrich(A, S, u, Fint, indexsetvec, res, params);
    } else {
        its = htawgm(A, S, u, Fint, indexsetvec, res, params2);
       // auto SF = applyScale(S, f, indexsetvec, 1e-06);
       // its = galerkin_pcg(A, S, u, f, indexsetvec, res,
       //                    true,
       //                    1e-08,
       //                    100,
       //                    1e-01,
       //                    1e-01,
       //                    1e-01,
       //                    1e-05);
    }

    std::cout << "Solver took " << its << " iterations to reach "
              << res << " accuracy" << std::endl;

    return 0;
}
