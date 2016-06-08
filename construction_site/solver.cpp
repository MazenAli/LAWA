#include <iostream>
#include <lawa/lawa.h>
#include <vector>
#include <cmath>

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
    Basis           basis(2);
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

    Integral1D cosi(basis, cosf, deltas, 4);
    Integral1D sini(basis, sinf, deltas, 4);
    Integral1D expi(basis, expf, deltas, 4);
    Integral1D x2i( basis, x2f,  deltas, 4);
    Integral1D x3i( basis, x3f,  deltas, 4);

    int rank = 2;
    int dim  = 4;
    SepCoeff        coeffs(rank, dim);
    SepCoeff        coeffset(rank, 1);
    IndexSetVec     indexsetvec(dim);
    lawa::SeparableFunctionD<T> F(fvec, rank, dim);
    MatInt                      derivs(rank, dim);
    for (int i=1; i<=rank; ++i) {
        for (int j=1; j<=dim; ++j) {
            derivs(i,j) = 0;//1;//(i+j)%2;
            _deltas.push_back(deltas);
        }
    }

    lawa::SeparableRHSD<T, Basis>   Fint(basis, F, _deltas, derivs);
    Index1D                         index1(0, 1, lawa::XWavelet);
    Index1D                         index2(2, 1, lawa::XWavelet);
    Index1D                         _index2(0, 3, lawa::XWavelet);
    Index1D                         index3(1, 2, lawa::XWavelet);
    Index1D                         index4(0, 1, lawa::XBSpline);
    Index1D                         index5(5, 3, lawa::XWavelet);
    std::vector<Index1D>            indexvec;
    std::vector<Index1D>            indexvec2;
    indexvec.push_back(index1);
    indexvec.push_back(index1);
    indexvec.push_back(index1);
    indexvec.push_back(index1);

    indexvec2.push_back(_index2);
    indexvec2.push_back(_index2);
    indexvec2.push_back(_index2);
    indexvec2.push_back(_index2);


    getFullIndexSet(basis, indexset, 0);
    indexset.insert(index3);
    getFullIndexSet(basis, indexset2, 1);

    std::cout << "The first index set is\n" << indexset
              << std::endl << "The second index set is\n" << indexset2
              << std::endl;

    IndexD                          indexmd(indexvec);
    IndexD                          indexmd2(indexvec2);

    indexsetvec[0] = indexset;
    for (int l=1; (unsigned)l<indexsetvec.size(); ++l) {
        indexsetvec[l] = indexset2;
    }

    genCoefficients(coeffs, Fint, indexsetvec);
    lawa::HTCoefficients<T, Basis>    f(dim, basis);
    lawa::HTCoefficients<T, Basis>    u(dim, basis);
    lawa::HTCoefficients<T, Basis>    r(dim, basis);
    set(f, coeffs);

    Laplace1D       LaplaceBil(basis);
    RefLaplace1D    RefLaplaceBil(basis.refinementbasis);
    Identity1D      IdentityBil(basis);
    RefIdentity1D   RefIdentityBil(basis.refinementbasis);
    LOp_Lapl1D      lapl(basis, basis, RefLaplaceBil, LaplaceBil);

    Sepop A(lapl, dim, dim);

    /* Testing eval by explicitly assembling matrix */
/*
    GeMat       Amat;
    DenseVector x;
    DenseVector Ax;
    unsigned max1d = maxintind(indexset, basis);
    unsigned size  = max1d*max1d*max1d*max1d;
    x.resize(size);
    Ax.resize(size);
    Amat.resize(size, size);

    for (auto& lambda1 : indexset) {
        for (auto& lambda2 : indexset) {
            for (auto& lambda3 : indexset) {
                for (auto& lambda4 : indexset) {
                    indexvec[0] = lambda1;
                    indexvec[1] = lambda2;
                    indexvec[2] = lambda3;
                    indexvec[3] = lambda4;
                    IndexD  lambdarow(indexvec);

                    unsigned k1 = maptoint(lambda1, basis);
                    unsigned k2 = maptoint(lambda2, basis);
                    unsigned k3 = maptoint(lambda3, basis);
                    unsigned k4 = maptoint(lambda4, basis);

                    unsigned row = (k1-1)*max1d*max1d*max1d
                                    +(k2-1)*max1d*max1d
                                    +(k3-1)*max1d
                                    +k4;

                    x(row) = f(lambdarow);

                    for (auto& lambdaj1 : indexset) {
                        for (auto& lambdaj2 : indexset) {
                            for (auto& lambdaj3 : indexset) {
                                for (auto& lambdaj4 : indexset) {

                                    unsigned j1 = maptoint(lambdaj1, basis);
                                    unsigned j2 = maptoint(lambdaj2, basis);
                                    unsigned j3 = maptoint(lambdaj3, basis);
                                    unsigned j4 = maptoint(lambdaj4, basis);

                                    unsigned col = (j1-1)*max1d*max1d*max1d
                                                    +(j2-1)*max1d*max1d
                                                    +(j3-1)*max1d
                                                    +j4;

                                    double a1, a2, a3, a4;

                                    a1 = LaplaceBil(lambda1, lambdaj1);

                                    if (k2==j2) a2=1.;
                                    else        a2=0.;

                                    if (k3==j3) a3=1.;
                                    else        a3=0.;

                                    if (k4==j4) a4=1.;
                                    else        a4=0.;

                                    Amat(row, col) = a1*a2*a3*a4;

                                    if (k1==j1) a1=1.;
                                    else        a1=0.;

                                    a2 = LaplaceBil(lambda2, lambdaj2);

                                    if (k3==j3) a3=1.;
                                    else        a3=0.;

                                    if (k4==j4) a4=1.;
                                    else        a4=0.;

                                    Amat(row, col) += a1*a2*a3*a4;

                                    if (k1==j1) a1=1.;
                                    else        a1=0.;

                                    if (k2==j2) a2=1.;
                                    else        a2=0.;

                                    a3 = LaplaceBil(lambda3, lambdaj3);

                                    if (k4==j4) a4=1.;
                                    else        a4=0.;

                                    Amat(row, col) += a1*a2*a3*a4;

                                    if (k1==j1) a1=1.;
                                    else        a1=0.;

                                    if (k2==j2) a2=1.;
                                    else        a2=0.;

                                    if (k3==j3) a3=1.;
                                    else        a3=0.;

                                    a4 = LaplaceBil(lambda4, lambdaj4);

                                    Amat(row, col) += a1*a2*a3*a4;
                                }
                            }
                        }
                    }

                }
            }
        }
    }

    Ax = Amat*x;


    std::cout << "A  =\n" << Amat  << std::endl;
    std::cout << "x  =\n" << x  << std::endl;
    std::cout << "Ax =\n" << Ax << std::endl;

    f.orthogonalize();
    u = eval(A, f, indexsetvec, indexsetvec, 1e-32);

    for (auto& lambda1 : indexset) {
        for (auto& lambda2 : indexset) {
            for (auto& lambda3 : indexset) {
                for (auto& lambda4 : indexset) {
                    indexvec[0] = lambda1;
                    indexvec[1] = lambda2;
                    indexvec[2] = lambda3;
                    indexvec[3] = lambda4;
                    IndexD  lambdarow(indexvec);

                    unsigned k1 = maptoint(lambda1, basis);
                    unsigned k2 = maptoint(lambda2, basis);
                    unsigned k3 = maptoint(lambda3, basis);
                    unsigned k4 = maptoint(lambda4, basis);

                    unsigned row = (k1-1)*max1d*max1d*max1d
                                    +(k2-1)*max1d*max1d
                                    +(k3-1)*max1d
                                    +k4;

                    x(row) = u(lambdarow);
                }
            }
        }
    }

    std::cout << "Ax by eval\n" << x << std::endl;

    Ax = Ax-x;
    std::cout << "Error = " << flens::blas::nrm2(Ax) << std::endl;

    exit(1);*/
    /* ---------------------------------------------- */

    lawa::Sepdiagscal<Basis>    S(dim, basis);
    double scaleps = 1e-02;
    setScaling(S, scaleps);

    std::cout << "***Testing solvers***\n\n";

    unsigned it;
    double residual;

    it = galerkin_pcg(A, S, u, f, indexsetvec, residual,
                      true, 1e-08, 5, 1e-01, 1e-09);

    std::cout << "\nScaling set last is\n" << S << std::endl;

    std::cout << "galerkin_pcg needed " << it << " iterations to reach r = "
              << residual << std::endl;

    exit(1);

    std::vector<IndexSet>    sweep(dim), total(dim);
    for (unsigned i=0; i<sweep.size(); ++i) {
        sweep[i]   = indexsetvec[i];
        total[i]   = indexsetvec[i];
    }

    std::cout << "+++Before eval total is--->\n";
    for (unsigned i=0; i<total.size(); ++i) {
        std::cout << "d = " << i+1 << std::endl;
        std::cout << total[i] << std::endl;
    }

    sweep = presidual(A, S, u, f, r, Fint, indexsetvec, sweep, total,
                         1e-10);

    T resnorm = nrm2(r);
    std::cout << "Approx residual is r = " << resnorm << std::endl;

    std::cout << "+++After eval total is--->\n";
    for (unsigned i=0; i<total.size(); ++i) {
        std::cout << "d = " << i+1 << std::endl;
        std::cout << total[i] << std::endl;
    }

    std::cout << "+++diff is--->\n";
    for (unsigned i=0; i<sweep.size(); ++i) {
        std::cout << "d = " << i+1 << std::endl;
        std::cout << sweep[i] << std::endl;
    }

    double alpha = 0.98;
    sweep = bulk(alpha, resnorm, residual, r, indexsetvec, sweep, total);
    std::cout << "alpha = " << alpha << std::endl;
    std::cout << "Actually the resulting ratio is "
              << nrm2(r)/resnorm << std::endl;

    std::cout << "+++new Lambda is--->\n";
    for (unsigned i=0; i<indexsetvec.size(); ++i) {
        std::cout << "d = " << i+1 << std::endl;
        std::cout << indexsetvec[i] << std::endl;
    }

    std::cout << "+++new sweep is--->\n";
    for (unsigned i=0; i<sweep.size(); ++i) {
        std::cout << "d = " << i+1 << std::endl;
        std::cout << sweep[i] << std::endl;
    }

    return 0;
}
