#include <iostream>
#include <lawa/lawa.h>


typedef     double                                  T;
typedef     flens::DenseVector<flens::Array<T>>     DenseVectorT;
typedef     flens::GeMatrix<flens::FullStorage
            <T, cxxblas::ColMajor>>                 FullColMatrixT;
typedef     lawa::Basis<T, lawa::Orthogonal,
            lawa::Interval, lawa::Multi>            Basis_XY;
typedef     lawa::TensorBasis2D<lawa::Adaptive,
            Basis_XY, Basis_XY>                     Basis2D;
typedef     lawa::SeparableRHS2D<T, Basis2D>        RHSINTEGRAL;
typedef     lawa::Index2D                           Index;
typedef     lawa::IndexSet<Index>                   IndexSet;
typedef     lawa::NoPreconditioner<T, Index>        NoPrec;
typedef     lawa::Coefficients<
            lawa::Lexicographical, T, Index>        DataType;
typedef     lawa::H1NormPreconditioner2D<T,
            Basis2D>                                Prec2D;


#define     a       1./35.
#define     b       1./35.
#define     mx      0.2
#define     my      0.2


double
xexp(double x)
{
    return std::exp(-1.*std::pow(x-mx, 2.)/(a*a));
}


double
yexp(double y)
{
    return std::exp(-1.*std::pow(y-my, 2.)/(b*b));
}


double
xexp_x(double x)
{
    return xexp(x)*(-1.*(x-mx)*2./(a*a));
}


int
main()
{
    // Reference RHS
    // Works only for functions with zero BCs
    DenseVectorT                    sing_support;
    FullColMatrixT                  nodeltas;
    lawa::SeparableFunction2D<T>    uref (xexp,   sing_support,
                                          xexp,   sing_support);
    lawa::SeparableFunction2D<T>    uxref(xexp_x, sing_support,
                                          xexp,   sing_support);
    lawa::SeparableFunction2D<T>    uyref(xexp,   sing_support,
                                          xexp_x, sing_support);


    int                             d(4);
    Basis_XY                        basis_xy(d, 2);
    Basis2D                         basis2d(basis_xy, basis_xy);
    RHSINTEGRAL                     _uh(basis2d, uref, nodeltas,
                                                       nodeltas, 100);
    NoPrec                          prec;
    lawa::RHS<T, Index,
    RHSINTEGRAL, NoPrec>            uh(_uh, prec);



    DataType Uh;
    T gamma = 0.2;
    T tol   = 1e-6;
    IndexSet Lambda;
    lawa::getSparseGridIndexSet(basis2d, Lambda, 1, 0, gamma);
    lawa::sample_f(basis2d, Lambda, uh, Uh, tol, true);
    T error = lawa::h1Error2D(basis2d, uref, uxref, uyref, Uh);
    std::cout << "Given L2 expected error " << tol
              << "\nH1 error is " << error << std::endl;
    return 0;
}
