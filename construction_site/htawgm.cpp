#include <iostream>
#include <lawa/lawa.h>
#include <vector>
#include <cmath>
#include <ctime>

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
x2(double x)
{
    return x*x;
}

double
x3(double x)
{
    return x*x*x;
}

double
mycos(double x)
{
    return std::cos(x);
}

    double
mysin(double x)
{
    return std::sin(x);
}

double
myexp(double x)
{
    return std::exp(x);
}


int
main()
{
    DenseVector     sings;
    GeMat           deltas(3,2);
    /*
    deltas(1,1) = 0.3;
    deltas(2,1) = 0.5;
    deltas(3,1) = 0.8;
    deltas(1,2) = 1.;
    deltas(2,2) = 2.;
    deltas(3,2) = 1.;
    */
    std::vector<GeMat> _deltas;
    Function        x2f(x2, sings);
    Function        x3f(x3, sings);
    Function        cosf(mycos, sings);
    Function        sinf(mysin, sings);
    Function        expf(myexp, sings);
    Basis           basis(4);
    basis.template enforceBoundaryCondition<lawa::DirichletBC>();
    IndexSet        indexset;
    IndexSet        indexset2;
    Coeff1D         coeff;

    std::vector<Function> fvec;
    fvec.push_back(cosf);
    fvec.push_back(sinf);
    fvec.push_back(sinf);
    fvec.push_back(expf);
    fvec.push_back(expf);
    fvec.push_back(x3f);
    fvec.push_back(x2f);
    fvec.push_back(expf);

    int rank = 2;
    int dim  = 4;
    SepCoeff        coeffs(rank, dim);
    SepCoeff        coeffs2(rank, dim);
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

    getFullIndexSet(basis, indexset, 4);
    //getFullIndexSet(basis, indexset2, 15);

    std::cout << "The index set size is\n" << indexset.size()
              << std::endl;
    //std::cout << "The index set size is\n" << indexset2.size()
    //          << std::endl;

    for (int l=0; (unsigned)l<indexsetvec.size(); ++l) {
        indexsetvec[l] = indexset;
        //indexsetvec2[l] = indexset2;
    }

    genCoefficients(coeffs, Fint, indexsetvec);
    genCoefficients(coeffs2, Fint, indexsetvec2);
    lawa::HTCoefficients<T, Basis>    f(dim, basis);
    lawa::HTCoefficients<T, Basis>    u(dim, basis);
    lawa::HTCoefficients<T, Basis>    r(dim, basis);
    set(f, coeffs);
/*
    set(u, coeffs2);

    double eps = 1e-06;
    DenseVector eps2(2);

    eps2(1) = eps*eps/2.;
    eps2(2) = eps*eps/2.;

    f.tree() = f.tree() + f.tree();
    f.tree() = f.tree() + f.tree();
    f.tree() = f.tree() + f.tree();
    f.tree() = f.tree() + f.tree();
    f.tree() = f.tree() + f.tree();
    f.tree() = f.tree() + f.tree();

    r = u;

    u.tree() = u.tree() + u.tree();
    u.tree() = u.tree() + u.tree();
    u.tree() = u.tree() + u.tree();
    u.tree() = u.tree() + u.tree();
    u.tree() = u.tree() + u.tree();
    u.tree() = u.tree() + u.tree();

    std::cout << "f rank = " << f.tree().max_rank() << std::endl;
    std::cout << "u rank = " << u.tree().max_rank() << std::endl;

    r = f;
    std::clock_t begin = std::clock();
    for (int k=1; k<=1; ++k) {
    f.truncate(eps);
    f = r;
    }
    std::clock_t end   = std::clock();

    std::cout << "Time to truncate f = "
              <<(double)(end-begin)/(double)CLOCKS_PER_SEC << std::endl;

    std::cout << "\n\n\n";

    r = u;
    begin = std::clock();
    for (int k=1; k<=1; ++k) {
    u.truncate(eps);
    u = r;
    }
    end   = std::clock();

    std::cout << "Time to truncate u = "
              <<(double)(end-begin)/(double)CLOCKS_PER_SEC << std::endl;

    exit(1);
*/
    Laplace1D       LaplaceBil(basis);
    RefLaplace1D    RefLaplaceBil(basis.refinementbasis);
    Identity1D      IdentityBil(basis);
    RefIdentity1D   RefIdentityBil(basis.refinementbasis);
    LOp_Lapl1D      lapl(basis, basis, RefLaplaceBil, LaplaceBil);

    Sepop A(lapl, dim, dim);

    lawa::Sepdiagscal<Basis>    S(dim, basis);
    setScaling(S, 0.5);
    S.set_nu(1e-01);

    lawa::HTAWGM_Params  params;
    params.maxit_pcg  = 200;
    params.maxit_awgm = 200;
    params.maxit_bulk = 10;
    params.tol_awgm   = 1e-08;
    params.trunc_pres = 1e-09;
    params.delta_pcg  = 1e-01;
    params.alpha      = 0.95;
    params.gamma      = 0.1;

    /*
    double cmax = 1./(1.+params.kappaA);
    cmax *= 0.9999;
    double est_acc   = cmax/(1.+cmax);
    double omega     = est_acc/(1.-est_acc);
    double alpha_max = (1.-omega)/std::sqrt(params.kappaA)-omega;
    double alpha     = (omega+alpha_max)/2.;*/

    /*
    double gamma = (1.-omega)*(alpha-omega)/((1.+omega)*params.kappaA);
    gamma *= 0.9999;*/

    std::cout << "HTAWGM params =\n";
    std::cout << params << std::endl;

    unsigned its;
    double   res;

/*
    its = galerkin_pcg2(A, S, u, f, indexsetvec, res,
                              true,
                              1e-08,
                              200,
                              1e-01,
                              1e-09);

    std::cout << "galerkin_pcg took " << its << " iterations to reach "
              << res << std::endl;
    std::cout << "Final S is set to\n" << S << std::endl;

    exit(1);*/

    its = htawgm(A, S, u, Fint, indexsetvec, res, params);

    std::cout << "htawgm took " << its << " iterations to reach "
              << res << " accuracy" << std::endl;
    std::cout << "Final scaling set to\n" << S << std::endl;

    return 0;
}
