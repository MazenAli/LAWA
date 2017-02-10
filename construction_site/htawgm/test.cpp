#include <iostream>
#include <lawa/lawa.h>
#include <lawa/lawa.h>
using namespace lawa;
using namespace std;
using namespace flens;

int
main()
{
    GeMatrix<FullStorage<double, ColMajor> > A(2, 2);
    A = 1, 1, 1, 1;

    flens::blas::scal(5., A(_, 1));

    cout << "A =\n" << A << endl;

    IndexSet<Index1D> set(50);
    index_eqfunction<Index1D> eq;
    Index1D lambda1(2, 3, XWavelet);
    Index1D lambda2(5, 3, XWavelet);

    set.insert(lambda1);
    set.insert(lambda2);

    cout << "Bucket count "         << set.bucket_count() << endl;
    cout << "Set size is "          << set.size()         << endl;
    cout << "Set load factor is "   << set.load_factor()  << endl;
    cout << "lambda1 is in bucket " << set.bucket(lambda1) << endl;
    cout << "lambda2 is in bucket " << set.bucket(lambda2) << endl;
    cout << "That bucket has size "
         << set.bucket_size(set.bucket(lambda1)) << endl;
    cout << "Hash value of lambda is " << set.hash_function()(lambda1) << endl;
    cout << "That bucket has size "
         << set.bucket_size(set.bucket(lambda2)) << endl;
    cout << "Hash value of lambda is " << set.hash_function()(lambda2) << endl;
    cout << "Accessing the bucket we get "
         << *set.begin(set.bucket(lambda1)) << endl;

    cout << "Are they equal? " << eq(*set.begin(set.bucket(lambda1)), lambda2)
                               << endl;

    return 0;
}
