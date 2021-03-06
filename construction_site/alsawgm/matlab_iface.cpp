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


int
main()
{
    int rank = 2;
    int dim  = 4;

    /* Generate example */
    Basis                       basis(2);
    basis.template              enforceBoundaryCondition<lawa::DirichletBC>();
    lawa::Mapwavind<Index1D>    map(dim);
    map.rehash(10);

    IndexSet                    indexset;
    IndexSet                    indexset2;
    IndexSet                    indexset3;
    getFullIndexSet(basis, indexset, 1);
    getFullIndexSet(basis, indexset2, 1);
    getFullIndexSet(basis, indexset3, 1);
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

    /* Test interface */
    rndinit(u, indexsetvec, 2);
    u.tree().orthogonalize();

    auto rhs = reduce_rhs(u, f);
    auto B   = reduce(A, u, indexsetvec, indexsetvec);
    std::vector<GeMat> x0(rhs.size()-1);
    IVector ranks(x0.size());
    extract_core(u.tree(), x0, ranks);

    std::cout << "Tree is\n";
    u.tree().print_w_UorB();
    std::cout << "\n";
    std::cout << "Cores are\n";
    for (unsigned k=0; k<x0.size(); ++k) {
        std::cout << "Core " << k+1 << "\n";
        std::cout << x0[k] << "\n";
        std::cout << "Rank : " << ranks(k+1) << "\n";
    }

    std::cout << "Initial matrices are\n";
    for (unsigned k=0; k<rhs.size(); ++k) {
        std::cout << "k = " << k+1 << std::endl;
        std::cout << "B = \n" << B[k] << std::endl;
        std::cout << "b = \n" << rhs[k] << std::endl;
    }

    /* Start MATLAB session */
    Engine *ep;
    std::cout << "Starting matlab...\n";
    if (!(ep = engOpen(""))) {
        std::cerr << "\nCan't start MATLAB engine\n" << std::endl;
        exit(1);
    }

    lawa::OptTTCoreParams params;
    std::vector<GeMat> output = lawa::optTTcoreLaplace(ep, B, rhs, x0, ranks, params);
    engClose(ep);
    insert_core(u.tree(), output, ranks);
    std::cout << "Result tree is\n";
    u.tree().print_w_UorB();

    return 0;
}
