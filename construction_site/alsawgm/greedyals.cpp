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
typedef     lawa::Sepop<
            lawa::AbstractLocalOperator1D<T,
            Basis>>
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


double
zero(double)
{
    return 0.;
}


int
main()
{
    int rank = 2;
    int dim  = 4;

    /* Generate example */
    Basis                       basis(2);
    basis.template              enforceBoundaryCondition<lawa::DirichletBC>();
    lawa::Mapwavind<Index1D>    map(dim);

    IndexSet                    indexset;
    IndexSet                    indexset2;
    IndexSet                    indexset3;
    getFullIndexSet(basis, indexset, 2);
    getFullIndexSet(basis, indexset2, 2);
    getFullIndexSet(basis, indexset3, 2);
    IndexSetVec     indexsetvec(dim);
    for (int l=0; (unsigned)l<indexsetvec.size(); ++l) {
        indexsetvec[l] = indexset;
    }
    indexsetvec[1] = indexset2;
    indexsetvec[3] = indexset3;

    double sp = 1.;
    HTCoeffTree                     f(dim, sp, basis, map);
    HTCoeffTree                     u(dim, sp, basis, map);
    SepCoeff                        coeffs(rank, dim);

    DenseVector             sings;
    Function                x2f(x2, sings);
    Function                x3f(x3, sings);
    Function                cosf(mycos, sings);
    Function                sinf(mysin, sings);
    Function                expf(myexp, sings);
    std::vector<Function>   fvec;
    fvec.push_back(cosf);
    fvec.push_back(sinf);
    fvec.push_back(sinf);
    fvec.push_back(expf);
    fvec.push_back(expf);
    fvec.push_back(x3f);
    fvec.push_back(x2f);
    fvec.push_back(expf);

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
    rndinit(u, indexsetvec, 1, 1e-02);

    for (int l=0; (unsigned)l<indexsetvec.size(); ++l) {
        std::cout << "Indexset " << l+1 << " : "
                  << indexsetvec[l].size() << std::endl;
    }

    lawa::Rank1UP_Params        p1;
    lawa::OptTTCoreParams       p2;
    lawa::GreedyALSParams       p3;
    lawa::H1NormPreconditioner1D<T, Basis> P(basis);
    //lawa::NoPreconditioner<T, Index1D>     P;
    double delta = 0.5;
    lawa::Sepdiagscal<Basis>    S(u.dim(), u.basis());
    setScaling(S, delta);
    S.set_nu(1e-02);

    /* Start MATLAB session */
    Engine *ep;
    if (!(ep = engOpen(""))) {
        std::cerr << "\nCan't start MATLAB engine\n" << std::endl;
        exit(1);
    }

    lawa::AgALSParams   params;
    params.maxit              = 10;
    params.gamma              = .5;
    params.rndinit            = 1e-02;
    params.r1update.sw        = true;
    params.r1update.orthog    = false;
    params.r1update.tol_als   = 1e-02;
    params.r1update.tol_cg    = 1e-08;
    params.r1update.check_res = false;
    params.greedyals.maxIt    = 25;

    std::cout << "The map is\n";
    for (const auto& mu : map.get_active()[0].left) {
        std::cout << mu.first << " mapsto " << mu.second << std::endl;
    }

    std::cout << "Solver parameters\n" << params << std::endl;
    double residual;
//    u.tree().orthogonalize();
//    (void) precrank1als_sym(A, P, u, f, indexsetvec, residual,
//                            true,
//                            false,
//                            true,
//                            1.,
//                            1e-03,
//                            10,
//                            1e-08,
//                            1e+02);
//    exit(1);
    unsigned it = agals_sym(ep, A, S, P, u, Fint, indexsetvec, residual, params);
//    params.coreopt.tol     = 1e-10;
//    params.coreopt.stag    = 1e-08;
//    params.greedyals.tol   = 1e-08;
//    params.r1update.sw     = true;
//    params.r1update.orthog = false;
//    unsigned it = greedyALS_sym(ep, A, P, u, f, indexsetvec, residual,
//                                params.r1update,
//                                params.coreopt,
//                                params.greedyals);
    std::cout << "AGALS took " << it << " iterations to reach relative residual "
              << residual << std::endl;

//    htucker::DimensionIndex idx(1);
//    idx[0] = 3;
//    auto toplot = extract(u, idx);
//    for (unsigned k=1; k<=toplot.rank(); ++k) {
//        std::cout << "Plotting function k=" << k << std::endl;
//        writeCoefficientsToFile(toplot(k, 1), k, "data_basisd3_nosworth");
//        std::string name = "basis_functiond3_nosworthlowcg_";
//        name            += std::to_string(k);
//        name            += ".dat";
//        plot(basis, toplot(k, 1), p, zero, zero, 0., 1.1, 1e-03, name.c_str());
//    }
//    std::cout << "Done...\n";

    engClose(ep);

    return 0;
}
