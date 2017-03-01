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

    rndinit(u, indexsetvec, 1);

    /* Some info */
    for (int j=1; j<=u.dim(); ++j) {
        htucker::DimensionIndex idx(1);
        idx[0]    = j;
        GeMat  Uj = extract(u.tree(), idx);
        std::cout << "Dimension             is " << j << "\n";
        std::cout << "The index set size    is " << indexsetvec[j-1].size() << std::endl;
        std::cout << "The linear     length is " << Uj.numRows()    << std::endl;
        std::cout << "The vectorized length is " << Uj.numRows()*Uj.numCols()
                  << "\n";
    }

    lawa::Rank1UP_Params    params;

    params.maxit_cg  = 100;
    params.max_sweep = 50;
    params.update    = false;
    params.tol_als   = 1e-02;
    double residual;
    unsigned nups    = 2;
    std::cout << "Update parameters are\n" << params;

    lawa::H1NormPreconditioner1D<T, Basis> P(basis);
    rank1greedy_sym(A, P, u, f, indexsetvec, residual, nups, params);

    auto rhs = reduce_rhs(u, f);
    auto B   = reduce(A, u, indexsetvec, indexsetvec);

    for (unsigned k=0; k<rhs.size(); ++k) {
        std::cout << "j = " << k+1 << std::endl;
        std::cout << "b = \n" << rhs[k] << std::endl;
        std::cout << "B = \n" << B[k] << std::endl;
    }

    return 0;
}
