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
    IndexSetVec     indexsetvec(dim*rank);
    IndexSetVec     indexsetvec2(dim);
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
    Index1D                         index3(1, 2, lawa::XWavelet);
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

    IndexSet    indexsetblub;
    getFullIndexSet(basis, indexsetblub, 1);
    IndexSet    indexset = indexsetblub; // testing copy
    indexset.insert(index3);
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

    indexsetvec2[0] = indexset;
    for (int l=1; (unsigned)l<indexsetvec2.size(); ++l) {
        indexsetvec2[l] = indexset2;
    }


    genCoefficients(coeffs, Fint, indexsetvec2);
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

    xpay(tree, idx, 1, -1., coeffset(2, 1));
    //axpy(tree, idx, 2, -1., coeffset(2, 1));
    //set(tree, idx, 1, 1.5, Fint(1, 1, indexset2));
    std::cout << "Post axpy tree is->\n";
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
    std::vector<lawa::AbstractLocalOperator1D<T, Basis>*>   ops_simple;
    for (int i=1; i<=dim; ++i) {
        ops_simple.push_back(&lapl);
        for (int j=1; j<=dim; ++j) {
            (j==i) ? ops.push_back(&lapl) : ops.push_back(&id);
        }
    }

    SepCoeff        cp(rank, dim);
    genCoefficients(cp, Fint, indexsetvec);
    lawa::HTCoefficients<T, Basis>    u(dim, basis);
    set(u, cp);
    u.orthogonalize();
    Sepop A(ops, dim, dim);
    Sepop As(ops_simple, dim, dim);
    Sepop Al(lapl, dim, dim);
    SepCoeff                          resultcp;
    SepCoeff                          resultcps;
    SepCoeff                          resultcpl;
    lawa::HTCoefficients<T, Basis>    resultu(dim, basis);
    lawa::HTCoefficients<T, Basis>    resultus(dim, basis);
    lawa::HTCoefficients<T, Basis>    resultul(dim, basis);
    lawa::HTCoefficients<T, Basis>    convertcp(dim, basis);
    lawa::HTCoefficients<T, Basis>    convertcps(dim, basis);
    lawa::HTCoefficients<T, Basis>    convertcpl(dim, basis);

    A.setrows(indexset, 1);
    A.setrows(indexset, 2);
    A.setrows(indexset, 3);
    A.setrows(indexset, 4);
    As.setrows(indexset, 1);
    As.setrows(indexset, 2);
    As.setrows(indexset, 3);
    As.setrows(indexset, 4);
    Al.setrows(indexset, 1);
    Al.setrows(indexset, 2);
    Al.setrows(indexset, 3);
    Al.setrows(indexset, 4);
    std::cout << "Return indexset is\n" << A.getrows(1) << std::endl;
    resultcp  = A*cp;
    resultcps = As*cp;
    resultcpl = Al*cp;
    resultu   = eval(A, u, 1e-32);
    resultus  = eval(As, u, 1e-32);
    resultul  = eval(Al, u, 1e-32);
    std::cout << "Result cp is ==>\n" << resultcp << std::endl;
    std::cout << "Result cps is ==>\n" << resultcps << std::endl;
    std::cout << "Result cpl is ==>\n" << resultcpl << std::endl;
    set(convertcp,  resultcp);
    set(convertcps, resultcps);
    set(convertcpl, resultcpl);

    std::cout << "Result cp in HT is ==>\n";
    convertcp.tree().print_w_UorB();
    std::cout << "Result cps in HT is ==>\n";
    convertcps.tree().print_w_UorB();
    std::cout << "Result cpl in HT is ==>\n";
    convertcpl.tree().print_w_UorB();

    std::cout << "\nResult u is ==>\n";
    resultu.tree().print_w_UorB();
    std::cout << "\nResult us is ==>\n";
    resultus.tree().print_w_UorB();
    std::cout << "\nResult ul is ==>\n";
    resultul.tree().print_w_UorB();

    htucker::HTuckerTree<T> diff = convertcp.tree()-resultu.tree();
    diff.orthogonalize();
    std::cout << "L2 error   = " << diff.L2normorthogonal() << std::endl;

    diff = convertcps.tree()-resultus.tree();
    diff.orthogonalize();
    std::cout << "L2 error s = " << diff.L2normorthogonal() << std::endl;

    diff = convertcpl.tree()-resultul.tree();
    diff.orthogonalize();
    std::cout << "L2 error l = " << diff.L2normorthogonal() << std::endl;

    diff = convertcp.tree()-resultul.tree();
    diff.orthogonalize();
    std::cout << "L2 error cp ul = "
              << diff.L2normorthogonal() << std::endl;

    diff = convertcps.tree()-resultul.tree();
    diff.orthogonalize();
    std::cout << "L2 error cps ul= "
              << diff.L2normorthogonal() << std::endl;

    diff = convertcp.tree()-resultus.tree();
    diff.orthogonalize();
    std::cout << "L2 error cp us= "
              << diff.L2normorthogonal() << std::endl;

    diff = resultu.tree()-resultus.tree();
    diff.orthogonalize();
    std::cout << "L2 error u us= "
              << diff.L2normorthogonal() << std::endl;

    diff = resultu.tree()-resultul.tree();
    diff.orthogonalize();
    std::cout << "L2 error u ul= "
              << diff.L2normorthogonal() << std::endl;

    diff = resultus.tree()-resultul.tree();
    diff.orthogonalize();
    std::cout << "L2 error us ul= "
              << diff.L2normorthogonal() << std::endl;

    resultu.tree().print_svs();

    std::cout << "\n\n***Checking scaling procedures\n\n";
    double order = 1.;
    set(u, cp);
    std::cout << "j0 = " << basis.j0 << std::endl;
    std::cout << "Omegamin2 for dim " << dim << " and order " << order
              << " is = " << compOmegamin2(basis, dim, order) << std::endl;
    std::cout << "The test vector in cp is " << cp << std::endl;
    lawa::Index1D testi(0, 3, lawa::XWavelet);
    std::cout << "The test index is " << testi << std::endl;
    std::cout << "Test values are "
              << std::setprecision(26) << cp(1, 3)[testi] << std::endl;
    std::cout << "and " << std::setprecision(26) << cp(2, 3)[testi]
              << std::endl;
    std::cout << "The indexscale for the cp vector is "
              << compIndexscale(cp, basis, order) << std::endl;
    std::cout << "The test vector in ht is\n";
    u.tree().print_w_UorB();
    std::cout << "The indexscale for the ht vector is "
              << compIndexscale(u, order) << std::endl;
    lawa::Sepdiagscal<Basis>    S(dim, basis);
    double scaleps = 1e-01;
    std::cout << "The scaling eps is set to " << scaleps << std::endl;
    setScaling(S, scaleps);
    std::cout << "The resulting params pre eval are\n" << S << std::endl;
    cp = S*cp;
    std::cout << "The resulting params post eval are\n" << S << std::endl;
    std::cout << "The result post scaling for cp is\n";
    std::cout << cp << std::endl;
    std::cout << "The results post scaling for ht is\n";
    //convertcp = S*convertcp;
    u = eval(S, u, 1e-08);
    std::cout << "The resulting params post 2nd eval are\n" << S << std::endl;
    u.tree().print_w_UorB();
    set(convertcpl,  cp);
    diff = u.tree()-convertcpl.tree();
    diff.orthogonalize();
    std::cout << "L2 error = " << diff.L2normorthogonal() << std::endl;

    std::cout << "Coeffs =>\n" << coeffs << std::endl;
    set(tree, coeffs);
    std::cout << "Tree =>\n";
    tree.tree().print_w_UorB();

    Index1D     indexnull(1, 1, lawa::XWavelet);
    Index1D     id2(1, 1, lawa::XWavelet);
    Index1D     id3(0, 1, lawa::XBSpline);
    Index1D     id4(1, 3, lawa::XWavelet);
    IndexD      idtest(dim);
    idtest(1) = indexnull;
    idtest(2) = id2;
    idtest(3) = id3;
    idtest(4) = id4;

    std::cout << "Test index is " << idtest << std::endl;
    std::cout << "First row number is " << maptoint(indexnull, basis)
              << std::endl;

    std::cout << "Coeffs evaled is " << coeffs(idtest) << std::endl;
    std::cout << "Tree evaled is " << tree(idtest) << std::endl;

    tree.orthogonalize();
    std::cout << "Tree evaled post orthog is " << tree(idtest) << std::endl;
    std::cout << "Tree post orthog is\n";
    tree.tree().print_w_UorB();

    std::cout << "\n\n---Testing new evals now---\n\n";

    IndexSetVec rows, cols;
    rows.push_back(indexset);
    rows.push_back(indexset2);
    rows.push_back(indexset2);
    rows.push_back(indexset2);
    cols.push_back(indexset);
    cols.push_back(indexset2);
    cols.push_back(indexset2);
    cols.push_back(indexset);

    cp.resize(rank, dim);
    genCoefficients(cp, Fint, indexsetvec2);
    std::cout << "The input cp is ==>\n" << cp << std::endl;

    int testrank = 2, testdim = 3;
    std::cout << "We will test application to rank "
              << testrank << " dim " << testdim << std::endl;
    Coeff1D     result;
    for (auto& row : rows[testdim-1]) {
        for (auto& col : cols[testdim-1]) {
            result[row] += LaplaceBil(row, col)*cp(testrank, testdim)[col];
        }
    }

    set(u, cp);
    u.orthogonalize();
    resultcp  = eval(A, cp, rows, cols);
    resultcps = eval(As, cp, rows, cols);
    resultcpl = eval(Al, cp, rows, cols);
    resultu   = eval(A, u, rows, cols, 1e-32);
    resultus  = eval(As, u, rows, cols, 1e-32);;
    resultul  = eval(Al, u, rows, cols, 1e-32);
    std::cout << "Result cp is ==>\n" << resultcp << std::endl;
    std::cout << "Result cps is ==>\n" << resultcps << std::endl;
    std::cout << "Result cpl is ==>\n" << resultcpl << std::endl;
    set(convertcp,  resultcp);
    set(convertcps, resultcps);
    set(convertcpl, resultcpl);

    Coeff1D delta = result - resultcp(6, 3);
    std::cout << "The result is\n" << result << std::endl;
    std::cout << "It should be\n" << resultcp(6, 3) << std::endl;
    std::cout << "The error is " << delta.norm() << std::endl;

    diff = convertcp.tree()-resultu.tree();
    diff.orthogonalize();
    std::cout << "L2 error cp u  = " << diff.L2normorthogonal() << std::endl;

    diff = convertcps.tree()-resultus.tree();
    diff.orthogonalize();
    std::cout << "L2 error cps us = " << diff.L2normorthogonal() << std::endl;

    diff = convertcpl.tree()-resultul.tree();
    diff.orthogonalize();
    std::cout << "L2 error cpl ul = " << diff.L2normorthogonal() << std::endl;

    diff = convertcps.tree()-resultu.tree();
    diff.orthogonalize();
    std::cout << "L2 error u cps   = " << diff.L2normorthogonal() << std::endl;

    diff = resultu.tree()-resultul.tree();
    diff.orthogonalize();
    std::cout << "L2 error u ul = " << diff.L2normorthogonal() << std::endl;

    diff = resultus.tree()-resultul.tree();
    diff.orthogonalize();
    std::cout << "L2 error us ul = " << diff.L2normorthogonal() << std::endl;

    diff = convertcp.tree()-convertcps.tree();
    diff.orthogonalize();
    std::cout << "L2 error cp cps = "
              << diff.L2normorthogonal() << std::endl;

    std::cout << "Now testing new scaling\n";
    resultcp  = eval(S, cp, indexsetvec2);
    resultcps = S*cp;
    set(convertcp, resultcp);
    set(convertcps, resultcps);

    diff = convertcps.tree()-convertcp.tree();
    diff.orthogonalize();
    std::cout << "Error 2 scalings = " << diff.L2normorthogonal() << std::endl;

    resultcp  = eval(S, cp, cols);
    set(convertcp, resultcp);
    std::cout << "Initial cp was " << cp << std::endl;
    std::cout << "Rescaled cp is" << resultcp << std::endl;

    resultu = eval(S, u, cols, 1e-32);
    diff = resultu.tree()-convertcp.tree();
    diff.orthogonalize();
    std::cout << "Does ht rescale do the same? -> "
              << diff.L2normorthogonal() << std::endl;

    std::cout << "The HT input u is\n";
    u.tree().print_w_UorB();
    std::cout <<"\n\nThe HT output is\n";
    resultu.tree().print_w_UorB();
    std::cout << "Interesting indexset " << indexset << std::endl;
    std::cout << "In integer indices these are==>\n";
    for (const auto& mu : indexset) {
        std::cout << "i = " << maptoint(mu, basis) << std::endl;
    }

    std::cout << "\n\nThe conversion is\n";
    convertcp.tree().print_w_UorB();
    std::cout << "\n\nThe conversion post orthog is\n";
    convertcp.orthogonalize();
    convertcp.tree().print_w_UorB();
    std::cout << "\n\nThe conversion post truncation is\n";
    convertcp.truncate(1e-32);;
    convertcp.tree().print_w_UorB();

    convertcp = resultu;
    resultu.truncate(1e-08);
    diff = convertcp.tree()-resultu.tree();
    diff.orthogonalize();
    std::cout << "Error truncate = " << diff.L2normorthogonal() << std::endl;

    resultu.truncate(3);
    std::cout << "Post first rank truncate the tensor looks like=>\n";
    resultu.tree().print_w_UorB();
    resultu.truncate(3);
    std::cout << "Post second rank truncate the tensor looks like=>\n";
    resultu.tree().print_w_UorB();
    resultu.truncate(3);
    std::cout << "Post third rank truncate the tensor looks like=>\n";
    resultu.tree().print_w_UorB();
    std::vector<DenseVector>    sigmas;
    resultu.orthogonalize_svd(sigmas);
    std::cout << "Post orthogonalize_svd the tensor looks like=>\n";
    resultu.tree().print_w_UorB();
    std::cout << "The singular values look like=>\n";
    for (unsigned int i=0; i<sigmas.size(); ++i) {
        std::cout << "Dim : " << i+1 << std::endl;
        std::cout << sigmas[i] << std::endl;
    }


    diff = convertcp.tree()-u.tree();
    diff.orthogonalize();
    std::cout << "Error truncate rank = " 
              << diff.L2normorthogonal() << std::endl;

    std::cout << "*\n*\n*\nTesting contractions ***\n";
    Coeff1D contractiond2 = contraction(resultu, indexset, sigmas[1], 2);
    std::cout << "The index set is\n" << indexset << std::endl;
    std::cout << "The contraction along dim 2 is " << contractiond2
              << std::endl;

    std::vector<Coeff1D> contractionall;
    contraction(resultu, indexsetvec2, sigmas, contractionall);
    for (unsigned int i=0; i<contractionall.size(); ++i) {
        std::cout << "Dim : " << i+1 << std::endl;
        std::cout << contractionall[i] << std::endl;
    }

    std::cout << "The ht tensor is\n";
    resultu.tree().print_w_UorB();
    SepCoeff    frame(3, dim);
    htucker::DimensionIndex idxx(1);
    idxx[0] = 1;
    frame = extract(resultu, idxx);
    std::cout << "The frame in dim 1 is\n" << frame << std::endl;

    for (const auto& mu : indexsetvec2[0]) {
        double entry = 0.;
        for (unsigned int j=1; j<=frame.rank(); ++j) {
            entry += sigmas[0](j)*frame(j, 1)[mu]*frame(j, 1)[mu];
        }
        std::cout << "For index " << mu << " the entry should be : "
                  << std::sqrt(entry) << std::endl;
    }

    return 0;
}


