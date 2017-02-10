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
one(double x)
{
    return 1.;
}

double
x2(double x)
{
    return 1.+x*x;
}

double
x3(double x)
{
    return 1.+x*x*x;
}

double
mycos(double x)
{
    return 1.+std::cos(x);
}

double
mysin(double x)
{
    return 1.+std::sin(x);
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

    std::vector<GeMat> _deltas;
    Function        onef(one, sings);
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
/*
    fvec.push_back(onef);
    fvec.push_back(onef);
    fvec.push_back(onef);
    fvec.push_back(onef);
*/
    int rank = 2;
    int dim  = 4;
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

    getFullIndexSet(basis, indexset, 2);
    Index1D mu1(2, 1, lawa::XWavelet);
    Index1D mu2(2, 2, lawa::XWavelet);
    Index1D mu3(16, 1, lawa::XWavelet);
    Index1D mu4(17, 8, lawa::XWavelet);
    Index1D mu5(17, 1, lawa::XWavelet);
    Index1D mu6(16, 8, lawa::XWavelet);
    Index1D mu7(4, 2, lawa::XWavelet);
    Index1D mu8(4, 11, lawa::XWavelet);
    Index1D mu9(5, 1, lawa::XWavelet);
    Index1D mu10(5, 2, lawa::XWavelet);
    Index1D mu11(5, 11, lawa::XWavelet);
    Index1D mu12(5, 21, lawa::XWavelet);
 /*   indexset.insert(mu1);
    indexset.insert(mu2);
    indexset.insert(mu3);
    indexset.insert(mu4);
    indexset.insert(mu5);
    indexset.insert(mu6);
    indexset.insert(mu7);
    indexset.insert(mu8);
    indexset.insert(mu9);
    indexset.insert(mu10);
    indexset.insert(mu11);
    indexset.insert(mu12);*/

    std::cout << "The index set size is\n" << indexset.size()
              << std::endl;

    for (int l=0; (unsigned)l<indexsetvec.size(); ++l) {
        indexsetvec[l] = indexset;
    }

    /* Map */
    lawa::Mapwavind<Index1D> map(dim);
    map.rehash(50);

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

    lawa::HTAWGM_Params  params;
    params.maxit_pcg  = 100;
    params.maxit_awgm = 100;
    params.tol_awgm   = 1e-08;
    params.delta1_pcg = 1e-01;
    params.delta2_pcg = 1e-01;
    params.delta3_pcg = 1e-01;
    params.alpha      = 0.95;
    params.recompr    = 1e-02;
    params.gamma      = 0.1;
    params.theta      = 1e-09;

    std::cout << "HTAWGM params =\n";
    std::cout << params << std::endl;

    unsigned its;
    double   res;

//    genCoefficients(coeffs, Fint, indexsetvec);
//    set(f, coeffs);
    setScaling(S, 0.5);

/*
    std::cout << "Scaling before\n";
    std::cout << S << std::endl;
    its = galerkin_pcg(A, S, u, f, indexsetvec, res, true, 1e-08,
                       100, 1e-02, 1e-02, 1e-02, 1e-09);
    std::cout << "Scaling after\n";
    std::cout << S << std::endl;
    std::cout << "It took " << its << " iterations to reach " << res
              << std::endl;
    exit(1);
*/
    its = htawgm(A, S, u, Fint, indexsetvec, res, params);

    std::cout << "htawgm took " << its << " iterations to reach "
              << res << " accuracy" << std::endl;
    std::cout << "Final scaling set to\n" << S << std::endl;

    return 0;
}