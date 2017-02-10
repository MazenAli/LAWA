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
    Index1D mu1(5, 3,  lawa::XWavelet);
    Index1D mu2(5, 5,  lawa::XWavelet);
    Index1D mu3(6, 10, lawa::XWavelet);
    Index1D mu4(6, 20, lawa::XWavelet);
    Index1D mu5(7, 8,  lawa::XWavelet);
    indexset.insert(mu1);
    indexset.insert(mu2);
    indexset.insert(mu3);
    indexset.insert(mu4);
    indexset.insert(mu5);
    //getFullIndexSet(basis, indexset2, 8);

    std::cout << "The index set size is\n" << indexset.size()
              << std::endl;
    //std::cout << "The index set size is\n" << indexset2.size()
    //          << std::endl;

    for (int l=0; (unsigned)l<indexsetvec.size(); ++l) {
        indexsetvec[l] = indexset;
        //indexsetvec2[l] = indexset2;
    }
    //indexsetvec[3] = indexset2;
    //indexsetvec2[3] = indexset2;

    /* Map */
    lawa::Mapwavind<Index1D> map(dim);
    map.rehash(500);
/*    map.rehash(3);
    Index1D  mu1(1, 1, lawa::XWavelet);
    Index1D  mu2(2, 2, lawa::XWavelet);
    Index1D  mu3(3, 3, lawa::XWavelet);
    Index1D  mu4(1, 2, lawa::XWavelet);
    Index1D  mu5(5, 3, lawa::XWavelet);
    std::cout << "buckets " << map.buckets(1) << std::endl;
    std::cout << "size " << map.size(1) << std::endl;

    std::cout << "index " << mu1 << " mapped to " << map(mu1, 1) << std::endl;
    std::cout << "sits in bucket " << map.get_active()[0].bucket(mu1) << std::endl;

    std::cout << "index " << mu2 << " mapped to " << map(mu2, 1) << std::endl;
    std::cout << "sits in bucket " << map.get_active()[0].bucket(mu2) << std::endl;

    std::cout << "index " << mu3 << " mapped to " << map(mu3, 1) << std::endl;
    std::cout << "sits in bucket " << map.get_active()[0].bucket(mu3) << std::endl;

    std::cout << "index " << mu4 << " mapped to " << map(mu4, 1) << std::endl;
    std::cout << "sits in bucket " << map.get_active()[0].bucket(mu4) << std::endl;

    std::cout << "index " << mu5 << " mapped to " << map(mu5, 1) << std::endl;
    std::cout << "sits in bucket " << map.get_active()[0].bucket(mu5) << std::endl;

    std::cout << "Mapping again...\n";

    std::cout << "index " << mu1 << " mapped to " << map(mu1, 1) << std::endl;
    std::cout << "sits in bucket " << map.get_active()[0].bucket(mu1) << std::endl;

    std::cout << "index " << mu2 << " mapped to " << map(mu2, 1) << std::endl;
    std::cout << "sits in bucket " << map.get_active()[0].bucket(mu2) << std::endl;

    std::cout << "index " << mu3 << " mapped to " << map(mu3, 1) << std::endl;
    std::cout << "sits in bucket " << map.get_active()[0].bucket(mu3) << std::endl;

    std::cout << "index " << mu4 << " mapped to " << map(mu4, 1) << std::endl;
    std::cout << "sits in bucket " << map.get_active()[0].bucket(mu4) << std::endl;

    std::cout << "index " << mu5 << " mapped to " << map(mu5, 1) << std::endl;
    std::cout << "sits in bucket " << map.get_active()[0].bucket(mu5) << std::endl;
    std::cout << "mapped back " << map(map(mu5, 1), 1) << std::endl;
*/

    genCoefficients(coeffs, Fint, indexsetvec);
    //genCoefficients(coeffs2, Fint, indexsetvec2);
    lawa::HTCoefficients<T, Basis>    f(dim, basis, map);
    lawa::HTCoefficients<T, Basis>    u(dim, basis, map);
    lawa::HTCoefficients<T, Basis>    r(dim, basis, map);
 /*   Index1D  mu1(1, 1, lawa::XWavelet);
    Index1D  mu2(1, 2, lawa::XWavelet);
    Index1D  mu3(3, 5, lawa::XWavelet);
    Index1D  mu4(4, 4, lawa::XWavelet);
    Index1D  mu5(5, 8, lawa::XWavelet);
    Index1D  mu6(6, 10, lawa::XWavelet);
    Index1D  mu7(7, 15, lawa::XWavelet);
    Index1D  mu8(8, 3, lawa::XWavelet);
    Index1D  mu9(9, 2, lawa::XWavelet);
    Index1D  mu10(10, 5, lawa::XWavelet);

    Coeff1D c1, c2, c3, c4, c5, c6, c7, c8;

    c1[mu1]  = 1.;
    c1[mu2]  = 1.;
    c1[mu3]  = 1.;

    c2[mu4]  = 1.;
    c2[mu5]  = 1.;

    c3[mu6]  = 1.;
    c3[mu7]  = 1.;
    c3[mu8]  = 1.;

    c4[mu9]  = 1.;
    c4[mu10] = 1.;

    c5[mu1] = 1.;
    c5[mu2] = 1.;
    c5[mu3] = 1.;

    c6[mu4]  = 1.;
    c6[mu5]  = 1.;

    c7[mu6]  = 1.;
    c7[mu7]  = 1.;
    c7[mu8]  = 1.;

    c8[mu9] = 1.;
    c8[mu10]  = 1.;

    indexsetvec[0].insert(mu1);
    indexsetvec[0].insert(mu2);
    indexsetvec[0].insert(mu3);

    indexsetvec[1].insert(mu4);
    indexsetvec[1].insert(mu5);

    indexsetvec[2].insert(mu6);
    indexsetvec[2].insert(mu7);
    indexsetvec[2].insert(mu8);

    indexsetvec[3].insert(mu9);
    indexsetvec[3].insert(mu10);

    coeffs(1, 1) = c1;
    coeffs(1, 2) = c2;
    coeffs(1, 3) = c3;
    coeffs(1, 4) = c4;

    coeffs(2, 1) = c5;
    coeffs(2, 2) = c6;
    coeffs(2, 3) = c7;
    coeffs(2, 4) = c8;*/

    set(f, coeffs);
    /*f.tree() = f.tree()+f.tree();
    double eps = 1e-08;
    f.tree().truncate_hsvd_soft(eps);
    u.tree() = f.tree();
    f.tree().truncate_hsvd_soft(eps);
    u.tree() = u.tree()-f.tree();
    u.tree().orthogonalize();
    f.tree().print_w_UorB();
    std::cout << std::endl;
    std::cout << "Soft truncation error = " << u.tree().L2normorthogonal() << std::endl;
    exit(1);*/
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
    S.set_nu(1e-04);

/*
    std::cout << "In CP form=>\n" << coeffs << std::endl;
    std::cout << "\n +++ In HT form +++\n";
    f.tree().print_w_UorB();
    std::cout << std::endl;

    u = eval(S, f, indexsetvec, 1e-08);
    std::cout << "Final scaling S set to\n" << S << std::endl;
    std::cout << "\n*\n*\nPost S eval\n*\n*\n";
    u.tree().print_w_UorB();

    u = evalS2(S, f, indexsetvec, 1e-08);
    std::cout << "Final scaling S2 set to\n" << S << std::endl;
    std::cout << "\n*\n*\nPost S2 eval\n*\n*\n";
    u.tree().print_w_UorB();

    exit(1);*/

    lawa::HTAWGM_Params  params;
    params.maxit_pcg  = 200;
    params.maxit_awgm = 200;
    params.tol_awgm   = 1e-08;
    params.delta1_pcg = 1e-02;
    params.delta2_pcg = 1e-02;
    params.delta3_pcg = 1e-01;
    params.alpha      = 0.9;
    params.recompr    = 0.25;
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

    its = galerkin_pcg(A, S, u, f, indexsetvec, res,
                              true,
                              1e-08,
                              100,
                              1e-01,
                              1e-01,
                              1e-01,
                              1e-09);

    std::cout << "galerkin_pcg took " << its << " iterations to reach "
              << res << std::endl;
    std::cout << "Final S is set to\n" << S << std::endl;

    exit(1);


    its = htawgm(A, S, u, Fint, indexsetvec, res, params);

    std::cout << "htawgm took " << its << " iterations to reach "
              << res << " accuracy" << std::endl;
    std::cout << "Final scaling set to\n" << S << std::endl;

    return 0;
}
