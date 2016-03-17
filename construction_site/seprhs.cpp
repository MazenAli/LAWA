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
    IndexSet        indexset, indexset2;
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
    IndexSetVec     indexsetvec(dim*rank);
    lawa::SeparableFunctionD<T> F(fvec, rank, dim);
    MatInt                      derivs(rank, dim);
    for (int i=1; i<=rank; ++i) {
        for (int j=1; j<=dim; ++j) {
            derivs(i,j) = 0;//1;//(i+j)%2;
            _deltas.push_back(deltas);
        }
    }
    std::cout << "derivs is " << derivs << std::endl;

    DenseVector x(4);
    x = 0., 0., 0., 0.;
    std::cout << "F at\n" << x << "\nis equal\n" << F.eval(x) << std::endl;
    x = 1., 1., 1., 1.;
    std::cout << "F at\n" << x << "\nis equal\n" << F.eval(x) << std::endl;
    x = 0.5, 1., 0.5, 1.;
    std::cout << "F at\n" << x << "\nis equal\n" << F.eval(x) << std::endl;
    x = 0.3, 0.2, 0.8, 0.7;
    std::cout << "F at\n" << x << "\nis equal\n" << F(x) << std::endl;
    x = 0.7, 1., 0., 0.33;
    std::cout << "F at\n" << x << "\nis equal\n" << F(x) << std::endl;

    lawa::SeparableRHSD<T, Basis>   Fint(basis, F, _deltas, derivs);
    Index1D                         index1(0, 1, lawa::XWavelet);
    Index1D                         index2(2, 1, lawa::XWavelet);
    Index1D                         _index2(0, 3, lawa::XWavelet);
    Index1D                         index3(3, 2, lawa::XWavelet);
    Index1D                         index4(0, 1, lawa::XBSpline);
    Index1D                         index5(5, 3, lawa::XWavelet);
    std::vector<Index1D>            indexvec;
    std::vector<Index1D>            indexvec2;
    indexvec.push_back(index1);
    indexvec.push_back(index1);
    indexvec.push_back(index1);
    indexvec.push_back(index1);
/*    indexvec.push_back(index2);
    indexvec.push_back(index3);
    indexvec.push_back(index4);*/

    indexvec2.push_back(_index2);
    indexvec2.push_back(_index2);
    indexvec2.push_back(_index2);
    indexvec2.push_back(_index2);

    /*
    indexset.insert(index1);
    indexset.insert(index2);
    indexset.insert(index3);
    indexset.insert(index4);

    indexset2.insert(index1);
    indexset2.insert(index2);
    indexset2.insert(index3);*/

    getFullIndexSet(basis, indexset, 1);
    getFullIndexSet(basis, indexset2, 2);

    IndexD                          indexmd(indexvec);
    IndexD                          indexmd2(indexvec2);

    std::cout << "Fint at\n" << indexmd << "\nis equal\n"
              << Fint.eval(indexmd) << std::endl;
    std::cout << "\nor\n" << Fint(indexmd) << std::endl;
    T check = Fint(indexmd);

    std::cout << "Test for lambda_1=" << index4 << std::endl;
    for (int i=1; i<=rank; ++i) {
        for (int j=1; j<=dim; ++j) {
            std::cout << "f_" << j << "^" << i << "="
                      << Fint(i, j, index4) << std::endl;
        }
    }

    int i = 2, j = 3;
    std::cout << "f_" << j << "^" << i << " should be "
              << x3i(index4) << std::endl;

    T res;
    res = cosi(index1)*sini(index2)*expi(index3)*x2i(index4)+
          sini(index1)*expi(index2)*x3i(index3)*expi(index4);
    std::cout << "Result should be " << res << std::endl;
    std::cout << "Error=" << res-check << std::endl;

    std::cout << "The indexset is\n" << indexset << std::endl;
    std::cout << "The coefficients are\n"
              << Fint(i, j, indexset) << std::endl;

    for (int l=0; (unsigned)l<indexsetvec.size(); ++l) {
        indexsetvec[l] = indexset;
    }

    genCoefficients(coeffs, Fint, indexsetvec);
    T check2 = coeffs(indexmd);
    std::cout << "The separable coefficients are\n" << coeffs << std::endl;

    std::cout << "Error sep coeff " << res-check2 << std::endl;

    std::cout << "2nd index is " << indexmd2 << std::endl;
    std::cout << "Coeffs is " << coeffs(indexmd2) << std::endl;
    std::cout << "But the rhs is " << Fint(indexmd2) << std::endl;

    lawa::HTCoefficients<T, Basis>    tree(dim, basis);

    std::cout << "\n\n***** HTCOEFFICIENTS *****\n\n";

    std::cout << "Evaluating empty tree -> does not work\n";

    set(tree, coeffs);

    std::cout << "Evaluating set tree at " << indexmd << " = \n"
              << tree(indexmd) << std::endl;
    std::cout << "Error is " << res-tree(indexmd) << std::endl;
    std::cout << "And at " << indexmd2 << " = \n"
              << tree(indexmd2) << std::endl;

    std::cout << "Checking maxintind for sets: \n" << indexset
              << "\nmax is\n"
              << lawa::maxintind(indexset, basis) << std::endl;

    std::cout << "\n***HT Tree -->**\n";
    tree.tree().print_w_UorB();

    htucker::DimensionIndex idx(1);
    idx[0] = 4;
    std::cout << "Index << " << idx << std::endl;
    std::cout << "idx min = " << idx.min() << ", idx max = " << idx.max()
              << std::endl;
    //indexset.insert(index5);
    set(tree, idx, 2, Fint(i, j, indexset));
    //std::cout << "Post set tree is->\n";
    //tree.tree().print_w_UorB();
    //std::cout << std::endl;

    coeffset(1, 1) = Fint(i, j, indexset);
    coeffset(2, 1) = Fint(1, 1, indexset2);

    lawa::HTCoefficients<T, Basis>    tree2(dim, basis);

    set(tree2, idx, coeffset);
    std::cout << "Post set mat tree is->\n";
    tree2.tree().print_w_UorB();

    //axpy(tree, idx, 1, -1., coeffset(2, 1));
    //axpy(tree, idx, 2, -1., coeffset(2, 1));
    //set(tree, idx, 1, 1.5, Fint(1, 1, indexset2));
    //std::cout << "Post axpy tree is->\n";
    //tree.tree().print_w_UorB();

    //axpy(tree, idx, 1., coeffset);
    //std::cout << "Post axpy mat tree is->\n";
    //tree.tree().print_w_UorB();

    std::cout << "Coeffset is " << coeffset << std::endl;

    std::cout << "Evaluating set tree at " << indexmd << " = \n"
              << tree(indexmd) << std::endl;

    std::cout << "Extracting 1st col... ->\n"
              << extract(tree2, idx, 1);
    std::cout << "Extracting all cols... ->\n"
              << extract(tree2, idx);


    /** Here come separable operators
     */


    std::cout << "\n/**\n  *\n  *\n  *Here are tests for separable operators"
              << "\n  *\n  *\n  *\n  */\n";
    Laplace1D       LaplaceBil(basis);
    RefLaplace1D    RefLaplaceBil(basis.refinementbasis);
    Identity1D      IdentityBil(basis);
    RefIdentity1D   RefIdentityBil(basis.refinementbasis);
    LOp_Lapl1D      lapl(basis, basis, RefLaplaceBil, LaplaceBil);
    LOp_Id1D        id(basis, basis, RefIdentityBil, IdentityBil);

    std::vector<lawa::AbstractLocalOperator1D<T, Basis>*>   ops;
    for (int i=1; i<=dim; ++i) {
        for (int j=1; j<=dim; ++j) {
            (j==i) ? ops.push_back(&lapl) : ops.push_back(&id);
        }
    }

    Sepop A(ops, dim, dim);
    lawa::HTCoefficients<T, Basis>    result(dim, basis);
    A.setIndexset(indexset);
    std::cout << "Return indexset is\n" << A.getIndexset() << std::endl;
    result = A*tree;
    std::cout << "Result is\n";
    result.tree().print_w_UorB();

    return 0;
}


