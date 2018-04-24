#include <iostream>
#include <flens/flens.h>
#include <lawa/lawa.h>

typedef double                                  T;
typedef flens::DenseVector<flens::Array<T> >    DenseVector;
typedef flens::DenseVector<flens::Array<int> >  IDV;

int
main()
{
    DenseVector x(4);
    x(1) = 0.1;
    x(2) = 1e-03;
    x(3) = 10.1;
    x(4) = 1.1;

    IDV idx;
    std::cout << "old x\n" << x << std::endl;
    flens::sort(x, idx);
    std::cout << "new x\n" << x << std::endl;
    std::cout << "idx\n" << idx << std::endl;

    return 0;
}

