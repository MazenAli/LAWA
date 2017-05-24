#include <iostream>
#include <iomanip>
#include <htucker/htucker.h>
#include <lawa/lawa.h>
#include <vector>
#include <cmath>
#include <cstdlib>

typedef double                                  T;
typedef htucker::HTuckerTree<T>                 HTTree;
typedef htucker::HTuckerTreeNode<T>             HNode;
typedef htucker::GeneralTreeIterator
        <HNode>                                 TreeIt;
typedef htucker::GeneralTreeNode<
        htucker::HTuckerTreeNode<T> >           GNode;
typedef lawa::Basis<T, lawa::Orthogonal,
        lawa::Interval, lawa::Multi>            Basis;
typedef lawa::Index1D                           Index1D;
typedef lawa::IndexSet<Index1D>                 IndexSet;
typedef std::vector<IndexSet>                   IndexSetVec;
typedef lawa::SepCoefficients<
        lawa::Lexicographical, T, Index1D>      SepCoeff;
typedef flens::GeMatrix<flens::FullStorage
        <T, cxxblas::ColMajor> >                GeMat;
typedef flens::GeMatrix<flens::FullStorage
        <int, cxxblas::ColMajor> >              MatInt;
typedef flens::DenseVector<flens::Array<T> >    DenseVector;
typedef lawa::Coefficients<
        lawa::Lexicographical, T,
        lawa::Index1D>                          Coeff1D;
typedef flens::DenseVector<
        flens::Array<unsigned> >                IVector;
typedef lawa::Function<T>                       Function;
typedef lawa::HTCoefficients<T, Basis>          HTCoeffTree;

typedef     lawa::AdaptiveLaplaceOperator1D<T,
            lawa::Orthogonal,
            lawa::Interval,
            lawa::Multi>                        Laplace1D;
typedef     lawa::AdaptiveLaplaceOperator1D<T,
            lawa::Orthogonal,
            lawa::Interval,
            lawa::MultiRefinement>              RefLaplace1D;
typedef     lawa::LocalOperator1D<Basis,
            Basis, RefLaplace1D, Laplace1D>     LOp_Lapl1D;
typedef     lawa::Sepop<LOp_Lapl1D>             Sepop;


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


double
zero(double)
{
    return 0.;
}


double
one(double)
{
    return 1.;
}


int
main(int argc, char* argv[])
{
    if (argc!=3) {
        std::cerr << "Usage: " << argv[0] << " dim level\n";
        return 1;
    }

    int dim  = atoi(argv[1]);
    int lev  = atoi(argv[2]);

    if (dim<=0) {
        std::cerr << "Dimension must be a positive integer\n";
        return 1;
    }
    int rank = 1;

    /* Generate example */
    Basis                       basis(2);
    basis.template              enforceBoundaryCondition<lawa::DirichletBC>();
    lawa::Mapwavind<Index1D>    map(dim);

    IndexSet                    indexset;
    IndexSet                    indexset2;
    IndexSet                    diff;
    getFullIndexSet(basis, indexset,  lev);
    getFullIndexSet(basis, indexset2,  lev-1);
    diff = indexset;
    for (auto& it : indexset2) {
        diff.erase(it);
    }
    IndexSetVec     indexsetvec(dim);
    IndexSetVec     indexsetvec2(dim);
    IndexSetVec     diffvec(dim);
    for (int l=0; (unsigned)l<indexsetvec.size(); ++l) {
        indexsetvec[l] = indexset;
        indexsetvec2[l] = indexset2;
        diffvec[l] = diff;
    }

    double sp = 1.;
    HTCoeffTree                     f(dim, sp, basis, map);
    HTCoeffTree                     u(dim, sp, basis, map);
    SepCoeff                        coeffs(rank, dim);

    DenseVector             sings;
    Function                onef(one, sings);
    std::vector<Function>   fvec;

    for (int i=0; i<dim; ++i) {
        fvec.push_back(onef);
    }

    lawa::SeparableFunctionD<T>     F(fvec, rank, dim);
    GeMat                           deltas(3,2);
    std::vector<GeMat>              _deltas;
    MatInt                          derivs(rank, dim);
    for (int i=1; i<=rank; ++i) {
        for (int j=1; j<=dim; ++j) {
            derivs(i,j) = 0;
            _deltas.push_back(deltas);
        }
    }
    lawa::SeparableRHSD<T, Basis>   Fint(basis, F, _deltas, derivs);
    genCoefficients(coeffs, Fint, indexsetvec);
    set(f, coeffs);
    RefLaplace1D    RefLaplaceBil(basis.refinementbasis);
    Laplace1D       LaplaceBil(basis);
    LOp_Lapl1D      lapl(basis, basis, RefLaplaceBil, LaplaceBil);
    Sepop           A(lapl, dim, dim);
    /* ---------------------------------------------------------- */

    /* Test greedy solver */
    rndinit(u, indexsetvec, 25, 1e-12);

    double a   = 0.7;
    double nrm = nrm2(u);
    auto   cp  = u;
    auto   cp2 = indexsetvec2;

    std::cout << "Best N\n";
    (void) bulkBestN(a, nrm, u, indexsetvec2, diffvec);
    std::cout << "Quasi best N\n";
    (void) bulk(a, nrm, cp, cp2, diffvec);

    exit(1);

    for (int l=0; (unsigned)l<indexsetvec.size(); ++l) {
        std::cout << "Indexset " << l+1 << " : "
                  << indexsetvec[l].size() << std::endl;
    }

    lawa::Rank1UP_Params                    p1;
    lawa::OptTTCoreParams                   p2;
    lawa::GreedyALSParams                   p3;
    lawa::DiagonalLevelPreconditioner1D<T>  P;
    lawa::NoPreconditioner<T, Index1D>      p;
    double delta = 0.5;
    lawa::Sepdiagscal<Basis>    S(u.dim(), u.basis());
    setScaling(S, delta);
    S.set_nu(1e-01);

    /* Start MATLAB session */
    Engine *ep;
    if (!(ep = engOpen("matlab -nojvm"))) {
        std::cerr << "\nCan't start MATLAB engine\n" << std::endl;
        exit(1);
    }

    lawa::AgALSParams   params;
    params.maxit              = 35;
    params.gamma              = 1e-01;
    params.r1update.update    = false;
    params.r1update.sw        = true;
    params.r1update.balance   = 500.;
    params.r1update.orthog    = true;
    params.r1update.tol_als   = 5e-02;
    params.r1update.tol_cg    = 1e-08;
    params.r1update.check_res = false;
    params.r1update.max_sweep = 20;
    params.r1update.verbose   = true;
    params.r1update.maxit_cg  = 500;
    params.greedyals.maxit    = 1;
    params.bulk               = 0.5;

    std::cout << "Solver parameters\n" << params << std::endl;
    double residual;
    auto start  = std::chrono::system_clock::now();
    unsigned it = agals_laplace(ep, A, S, u, Fint, indexsetvec, residual, params);
//    params.r1update.sw        = true;
//    params.r1update.balance   = 500.;
//    params.r1update.orthog    = true;
//    params.r1update.tol_als   = 1e-02;
//    params.r1update.max_sweep = 50;
//    params.r1update.maxit_cg  = 150;
//    params.coreopt.tol        = 1e-08;
//    params.coreopt.stag       = 1e-07;
//    params.greedyals.tol      = 1e-07;
//    params.greedyals.stag     = 1e-07;
//    params.greedyals.maxit    = 30;
//    auto start  = std::chrono::system_clock::now();
//    unsigned it = greedyALS_laplace(ep, A, u, f, indexsetvec, residual,
//                                    params.r1update,
//                                    params.coreopt,
//                                    params.greedyals);
    auto end     = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start);
    std::cout << "AGALS took " << it << " iterations to reach relative residual "
            << residual << std::endl;
    std::cout << "It took " << elapsed.count() << " secs\n";

    htucker::DimensionIndex idx(1);
    idx[0] = 1;
    auto toplot = extract(u, idx);
    for (unsigned k=1; k<=toplot.rank(); ++k) {
        std::cout << "Plotting function k=" << k << std::endl;
        //writeCoefficientsToFile(toplot(k, 1), k, "data_basisd1");
        std::string name = "basis_functiond1_";
        name            += std::to_string(k);
        name            += ".dat";
        plot(basis, toplot(k, 1), p, zero, zero, 0., 1.1, 1e-03, name.c_str());
    }
    std::cout << "Done...\n";

    engClose(ep);

    return 0;
}
