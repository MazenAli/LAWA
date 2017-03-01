#include <iostream>
#include <iomanip>
#include <htucker/htucker.h>
#include <lawa/lawa.h>
#include <vector>
#include <cmath>

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

#include <engine.h>
#include <cstring>
std::vector<std::vector<GeMat>>
interface(const std::vector<std::vector<GeMat>>& B);


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
    getFullIndexSet(basis, indexset, 2);
    getFullIndexSet(basis, indexset2, 2);
    getFullIndexSet(basis, indexset3, 1);
    IndexSetVec     indexsetvec(dim);
    for (int l=0; (unsigned)l<indexsetvec.size(); ++l) {
        indexsetvec[l] = indexset;
    }
    indexsetvec[1] = indexset2;
    indexsetvec[3] = indexset3;


    HTCoeffTree                     f(dim, basis, map);
    HTCoeffTree                     u(dim, basis, map);
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
    rndinit(u, indexsetvec, 3);

    auto rhs = reduce_rhs(u, f);
    auto B   = reduce(A, u, indexsetvec, indexsetvec);

    std::vector<std::vector<GeMat>> input(2);
    std::vector<std::vector<GeMat>> output(2);
    input[0].resize(rhs.size());
    input[1].resize(rhs.size());


    std::cout << "Initial matrices are\n";
    for (unsigned k=0; k<rhs.size(); ++k) {
        std::cout << "k = " << k+1 << std::endl;
        std::cout << "B = \n" << B[k] << std::endl;
        std::cout << "b = \n" << rhs[k] << std::endl;
        input[0][k] = B[k];
        input[1][k] = rhs[k];
    }

    /* Call here */
    /* Start MATLAB session */
    Engine *ep;
    std::cout << "Starting matlab...\n";
    if (!(ep = engOpen(""))) {
        std::cerr << "\nCan't start MATLAB engine\n" << std::endl;
        exit(1);
    }

    output = lawa::optcore(ep, input);
    engClose(ep);
    std::cout << "Return matrices are\n";
    for (unsigned k=0; k<output.size(); ++k) {
        for (unsigned j=0; j<output[k].size(); ++j) {
            std::cout << "k = " << k+1 << ", j= " << j+1 << "\n";
            std::cout << output[k][j] << "\n";
        }
    }

    return 0;
}


std::vector<std::vector<GeMat>>
interface(const std::vector<std::vector<GeMat>>& B)
{
    mxArray *icell = flens::matlab::create2DCell(B);

    /* Start MATLAB session */
    Engine *ep;
    std::cout << "Starting matlab...\n";
    if (!(ep = engOpen(""))) {
        std::cerr << "\nCan't start MATLAB engine\n" << std::endl;
        exit(1);
    }

    /* Pass cell */
    engPutVariable(ep, "RedSys", icell);

    /* Call function */
    std::cout << "Displaying contents from MATLAB\n";
    engEvalString(ep, "disp(RedSys)");
    engEvalString(ep, "X = RedSys");

    /* Process result (change later) */
    mxArray *ocell = nullptr;
    ocell          = engGetVariable(ep, "X");
    if (!ocell) {
        std::cerr << "interface: Reading variable from workspace failed\n";
    }

    std::vector<std::vector<GeMat>> ret = flens::matlab::read2DCell(ocell);

    /* Clean up, end session */
    mxDestroyArray(icell);
    mxDestroyArray(ocell);
    engClose(ep);

    return ret;
}
