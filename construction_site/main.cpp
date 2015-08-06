#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <flens/flens.cxx>
#include <htucker/htucker.h>
#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/constructions/basis.h>
#include <lawa/constructions/interval/multi/basis.h>
#include <lawa/methods/adaptive/operators/pdeoperators1d/adaptivehelmholtzoperatoroptimized1d.h>
#include <construction_site/coeffframe.h>
#include <construction_site/htcoeffnode.h>

typedef double                                                      T;
typedef lawa::Basis<double, lawa::Orthogonal, lawa::Interval,
            lawa::Multi>                                            MBasis;
typedef lawa::HTCoeffNode<T, MBasis, lawa::Index1D>                 HTCoeffNode;
typedef lawa::CoeffFrame<lawa::Lexicographical, T, lawa::Index1D>   HTCoeffFrame;
typedef lawa::AdaptiveHelmholtzOperatorOptimized1D<T,
        lawa::Orthogonal, lawa::Interval, lawa::Multi>              HeHoOp;
typedef lawa::Coefficients<lawa::Lexicographical,T,lawa::Index1D>   Coefficients;

MBasis   basis(4);
htucker::DimensionIndex di(3, 1);


void
apply(HTCoeffNode& u, const T eps)
{
    T c = 0;
    HeHoOp op(basis, c);

    int numCols = u.numCols();
    HTCoeffFrame U(numCols), ret(numCols);

    typedef typename HTCoeffNode::const_iterator const_it;

    for (int j=1; j<=numCols; ++j) {
        for (const_it it=u.cbegin(j); it!=u.cend(j); ++it) {
            U((*it), j) = u((*it), j);
        }
        op.apply(U[j], eps, ret[j]);
    }

    u.setCoeff(ret);
}


void
addCoefficients(Coefficients& w_p, const int p,
                lawa::Coefficients<lawa::Bucket,T,lawa::Index1D>& v_bucket,
                const HTCoeffNode& U, const int col)
{
     typedef typename   lawa::Coefficients<lawa::Bucket,T,
                        lawa::Index1D>::BucketEntry::const_iterator const_it;

     for (const_it  it=v_bucket.buckets[p].begin();
                    it!=v_bucket.buckets[p].end(); ++it) {
        if (U.isActive((**it).first, col)) {
            w_p[(**it).first] = U((**it).first, col);
        }
    }
}


void
compContrac(const int i, Coefficients& pi,
            /*to be replaced integrated into class*/const HTCoeffNode& U)
{
    typedef typename HTCoeffNode::const_iterator    active_const_it;
    typedef typename Coefficients::iterator         coeff_it;

    // Use i to determine node...
    int dummy = i+1;

    // TODO: Check here if sigmas are squared already!
    for (int i=1; i<=U.numCols(); ++i) {
        T sigma = U.getNode().getSigma()(i);
        for (active_const_it it=U.cbegin(i); it!=U.cend(i); ++it) {
            pi[*it] = std::pow(U(*it, i)*sigma, 2.);
        }
    }

    for (coeff_it it=pi.begin(); it!=pi.end(); ++it) {
        pi[(*it).first] = std::sqrt((*it).second);
    }
}


void
apply_i(HTCoeffNode& U, const T eps, const Coefficients& v)
{
    typedef typename Coefficients::const_iterator       const_coeff1d_it;
    typedef typename Coefficients::iterator             coeff1d_it;
    typedef typename lawa::IndexSet<lawa::Index1D>::const_iterator
                                                        const_set1d_it;

    T c = 0.;
    HeHoOp op(basis, c);
    T comprConst = op.CA;
    T convFactor = std::min(basis.d-1.5,1.5);

    lawa::Coefficients<lawa::Bucket,T,lawa::Index1D> v_bucket;
    T tol = 0.5*eps/comprConst;
    v_bucket.bucketsort(v,tol);

    long double squared_v_norm = (long double)std::pow(v.norm(2.),2.);
    long double squared_v_bucket_norm = 0.;
    T delta=0.;

    int l=0;
    int support_size_all_buckets=0;
    int numCols = U.numCols();
    HTCoeffFrame ret(numCols);

    for (int i=0; i<(int)v_bucket.buckets.size(); ++i) {
        squared_v_bucket_norm +=    v_bucket.bucket_ell2norms[i]*
                                    v_bucket.bucket_ell2norms[i];
        T squared_delta = fabs(squared_v_norm - squared_v_bucket_norm);
        support_size_all_buckets += v_bucket.buckets[i].size();
        delta = std::sqrt(squared_delta);
        ++l;
        if (squared_delta<tol*tol) {
            break;
        }
    }

    if (delta>eps) delta = eps/2.; // round off errors

    for (int k=1; k<=numCols; ++k) {
        printf("\n###Iterating over column %d###\n", k);
        for (int i=0; i<l; ++i) {
            Coefficients w_p;
            addCoefficients(w_p, i, v_bucket, U, k);
            if (w_p.size()==0) continue;

            T numerator = w_p.norm(2.) * support_size_all_buckets;
            T denominator = w_p.size() * (eps-delta) / comprConst;
            int jp =    (int)std::max(
                        (std::log(numerator/denominator) /
                        std::log(2.)) / convFactor, 0.);

            printf("\n###Iterating over bucket %d###\n", i);
            for (const_coeff1d_it it=w_p.begin(); it!=w_p.end(); ++it) {
                lawa::Index1D colindex = (*it).first;
                T prec_colindex = op.prec(colindex);
                lawa::IndexSet<lawa::Index1D> Lambda_v;
                Lambda_v =  lambdaTilde1d_PDE(colindex, basis,jp, basis.j0,
                            std::min(colindex.j+jp,25), false);
                std::cout   << "\n###Iterating over colindex "
                            << colindex << " ###\n";
                for (const_set1d_it mu=Lambda_v.begin();
                                    mu!=Lambda_v.end(); ++mu) {
                    std::cout   << "###Active column entry is "
                                << *mu << " ###\n";
                    ret(*mu, k) +=  op.laplace_data1d(*mu, colindex)*
                                    prec_colindex * (*it).second;
                }
            }
        }

        for (coeff1d_it it=ret[k].begin(); it!=ret[k].end(); ++it) {
            (*it).second *= op.prec((*it).first);
        }
    }

    U.setCoeff(ret);
}


void
buildCoeff(Coefficients& v, const HTCoeffNode& U, const int col)
{
    typedef typename HTCoeffNode::const_iterator const_it;

    for (const_it it=U.cbegin(col); it!=U.cend(col); ++it) {
        v[(*it)] = U((*it), col);
    }
}


void
testApply(int j, int rank)
{
    srand(0);
    htucker::HTuckerTreeNode<T> htnode(di), htcopy(di), httemp(di);
    HTCoeffNode ht(htnode, basis), copy(htnode, basis), temp(htnode, basis);
    HTCoeffFrame u(rank+1);
    printf("\n\n---------Testing apply upto j=%d----------------\n", j);
    printf("---------           rank=%d             ----------------\n\n", rank);

    int j0 = basis.j0;
    bool flag(false);
    int offset(0);
    for (int l=1; l<=rank; ++l) {
        for (int k=basis.mra.rangeI(j0).firstIndex();
                 k<=basis.mra.rangeI(j0).lastIndex(); ++k) {
            offset =    basis.mra.rangeI(j0).lastIndex()-2*k
                        +basis.mra.rangeI(j0).firstIndex();
            lawa::XType type = lawa::XBSpline;
            if (l==1) offset = 0;
            lawa::Index1D lambda(j0, k+offset, type);
            double r = ((double) rand() / (double) RAND_MAX);
            u(lambda, l) = r;
        }
        for (int i=basis.j0; i<=j; ++i) {
            for (int k=basis.rangeJ(i).firstIndex();
                     k<=basis.rangeJ(i).lastIndex(); ++k) {
                offset =    basis.rangeJ(i).lastIndex()-2*k
                            +basis.rangeJ(i).firstIndex();
                lawa::XType type = lawa::XWavelet;
                if (l==1) offset = 0;
                lawa::Index1D lambda(i, k+offset, type);
                double r = ((double) rand() / (double) RAND_MAX);
                u(lambda, l) = r;
                /*if (l!=2 || !flag) {
                    u(lambda, l) = count/scale;
                    flag = true;
                } else {
                    u(lambda, l) = 0.1*(double)count;
                    flag = false;
                }*/
            }
        }
    }

    ht.setCoeff(u);
    copy.setCoeff(u);
    temp.setCoeff(u);

    printf("Before apply()\n");
    std::cout << ht.getFrame() << std::endl;
    T eps = .1;
    apply(ht, eps);
    printf("After apply()\n");
    std::cout << ht.getFrame() << std::endl;

    printf("\nNow comparing with apply_i\n\n");
    for (int i=1; i<=copy.numCols(); ++i) {
        printf("\n\n\n---###---Using column %d as reference\n---###---", i);
        Coefficients v;
        temp.setCoeff(u);
        buildCoeff(v, temp, i);
        apply_i(temp, eps, v);
        printf("\n\n\n---###---Using column %d as reference\n---###---", i);
        std::cout << temp.getFrame() << std::endl;
    }
}


int
main()
{
    htucker::HTuckerTreeNode<T> ht(di);
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  A(3,3);
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  B(4,4);
    flens::DenseVector<flens::Array<T> >                        v1(2);
    flens::DenseVector<flens::Array<T> >                        v2(3);
    ht.setUorB(A);
    ht.setUorB(B);
    ht.setSigma(v1);
    ht.setSigma(v2);
    std::cout << ht.getUorB() << std::endl;
    std::cout << ht.getSigma() << std::endl;

    testApply(3, 5);
    std::cerr << "This is a sample error output\n";
    return 0;
}
