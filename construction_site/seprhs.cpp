#include <iostream>
#include <lawa/lawa.h>
#include <vector>
#include <cmath>

typedef double                                  T;
typedef flens::DenseVector<flens::Array<T> >    DenseVector;
typedef lawa::Function<T>                       Function;

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
    Function        x2f(x2, sings);
    Function        x3f(x3, sings);
    Function        cosf(mycos, sings);
    Function        sinf(mysin, sings);
    Function        expf(myexp, sings);

    std::vector<Function> fvec;
    fvec.push_back(cosf);
    fvec.push_back(sinf);
    fvec.push_back(sinf);
    fvec.push_back(expf);
    fvec.push_back(expf);
    fvec.push_back(x3f);
    fvec.push_back(x2f);
    fvec.push_back(expf);

    int rank = 2;
    int dim  = 4;
    lawa::SeparableFunctionD<T> F(fvec, rank, dim);

    DenseVector x(4);
    x = 0., 0., 0., 0.;
    std::cout << "F at\n" << x << "\nis equal\n" << F(x) << std::endl;
    x = 1., 1., 1., 1.;
    std::cout << "F at\n" << x << "\nis equal\n" << F(x) << std::endl;
    x = 0.5, 1., 0.5, 1.;
    std::cout << "F at\n" << x << "\nis equal\n" << F(x) << std::endl;
    x = 0.3, 0.2, 0.8, 0.7;
    std::cout << "F at\n" << x << "\nis equal\n" << F(x) << std::endl;
    x = 0.7, 1., 0., 0.33;
    std::cout << "F at\n" << x << "\nis equal\n" << F(x) << std::endl;

    return 0;
}


