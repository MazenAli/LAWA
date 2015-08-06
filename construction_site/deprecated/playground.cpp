#include <iostream>
#include <flens/flens.cxx>
#include <htucker/htucker.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>

typedef double                                                      T;
typedef flens::GeMatrix<flens::FullStorage<T, flens::ColMajor> >    GeMatrix;
typedef flens::DenseVector<flens::Array<T> >                        DenseVector;
typedef htucker::HTuckerTreeNode<T>                                 HTNode;
typedef lawa::Coefficients< lawa::Lexicographical,
                            double, lawa::Index1D>                  Coefficients;
const   flens::Underscore<GeMatrix::IndexType>                      _;


int
main()
{
    GeMatrix            A(3, 3), C(5,3);
    DenseVector::View   v = A(_, 1);
    htucker::DimensionIndex j(1);
    j.setValue(1, 3);
    HTNode              U(j);

    Coefficients        u;
    lawa::Index1D       l1(0, 2, lawa::XBSpline),
                        l2(1, 3, lawa::XWavelet),
                        l3(4, 1, lawa::XWavelet),
                        l4(4, 5, lawa::XWavelet),
                        l5(0, -1, lawa::XBSpline),
                        l6(8, -1, lawa::XWavelet);
    u[l1] = -1.;
    u[l2] = 0.5;
    u[l3] = 1.7;
    u[l4] = -2.1;
    u[l5] = 2.3;
    u[l6] = 0.3;

    printf("Printing coefficient values...\n");
    std::cout << u << std::endl;

    A = 1, 2, 3,
        4, 5, 6,
        7, 8, 9;
    U.setUorB(A);
    C(_(1,3),_) = A;
    GeMatrix::ConstView B = U.getUorB();

    std::cout << "B = " << B << std::endl;
    std::cout << "A = " << A << std::endl;
    std::cout << "C = " << C << std::endl;

    return 0;
}
