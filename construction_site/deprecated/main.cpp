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
printBasis(const MBasis& _basis, int j)
{
    printf("Printing Basis upto level j = %d\n\n", j);
    int j0 = _basis.j0;
    int mracard, card;
    int kl, kr, mkl, mkr;

    for (int i=j0; i<=j; ++i) {
        mracard = _basis.mra.cardI(i);
        card    = _basis.cardJ(i);
        kl      = _basis.rangeJ(i).firstIndex();
        kr      = _basis.rangeJ(i).lastIndex();
        mkl     = _basis.mra.rangeI(i).firstIndex();
        mkr     = _basis.mra.rangeI(i).lastIndex();

        printf("On level j = %d, card = %d, mracard = %d\n", i, card, mracard);
        printf("The range is between %d - %d\n", kl, kr);
        printf("The MRA range is between %d - %d\n", mkl, mkr);
    }
}

void
testMap(int j)
{
    htucker::HTuckerTreeNode<T> htnode(di);
    HTCoeffNode ht(htnode, basis);
    printf("\n\n---------Testing Map upto j=%d----------------\n\n", j);

    int j0 = basis.j0;
    for (int k=basis.mra.rangeI(j0).firstIndex();
             k<=basis.mra.rangeI(j0).lastIndex(); ++k) {
        lawa::XType type = lawa::XBSpline;
        lawa::Index1D lambda(j0, k, type);
        int coeff = lawa::mapCoeff(lambda, basis);
        printf("Lambda is (j,k,xtype) = (%d,%d,%d)\n", j0, k, (int)type);
        printf("Mapped number is coeff = %d\n", coeff);
    }

    printf("\n-------------\n");
    for (int i=basis.j0; i<=j; ++i) {
        for (int k=basis.rangeJ(i).firstIndex();
                 k<=basis.rangeJ(i).lastIndex(); ++k) {
            lawa::XType type = lawa::XWavelet;
            lawa::Index1D lambda(i, k, type);
            int coeff = lawa::mapCoeff(lambda, basis);
            printf("Lambda is (j,k,xtype) = (%d,%d,%d)\n", i, k, (int)type);
            printf("Mapped number is coeff = %d\n", coeff);
        }
    }
}


void
testIterators(HTCoeffNode& ht)
{
    printf("\n\n ------------------- Testing Iterators -------------------\n\n");

    typedef typename HTCoeffNode::const_iterator    const_it;
    typedef typename HTCoeffNode::iterator          iter;
    int numCols = ht.getRank();

    for (int j=1; j<= numCols; ++j) {
        printf("Column j=%d\n", j);
        for (iter it=ht.begin(j); it!=ht.end(j); ++it) {
            std::cout   << "u(["<< (*it) << "], " << j
                        << ") = " << ht((*it), j) << std::endl;
        }
    }

    printf("\n -------------- DONE --------------\n");
}


void
testHTCoeffNode(int j, int rank)
{
    htucker::HTuckerTreeNode<T> htnode(di);
    HTCoeffNode ht(htnode, basis);
    HTCoeffFrame u(rank);
    printf("\n\n---------Testing HTCoeffNode upto j=%d----------------\n", j);
    printf("---------           rank=%d             ----------------\n\n", rank);

    int count(0);
    int j0 = basis.j0;
    for (int l=1; l<=rank; ++l) {
        for (int k=basis.mra.rangeI(j0).firstIndex();
                 k<=basis.mra.rangeI(j0).lastIndex(); ++k) {
            lawa::XType type = lawa::XBSpline;
            lawa::Index1D lambda(j0, k, type);
            ++count;
            u(lambda, l) = count;
        }
        printf("Setting coefficients for scaling fncs\n");
        ht.setCoeff(u);
        std::cout << ht.getFrame() << std::endl;

        printf("\n...Printing Active Coefficients...\n");
        ht.printActive();

        printf("\n-------------\n");
        for (int i=basis.j0; i<=j; ++i) {
            for (int k=basis.rangeJ(i).firstIndex();
                     k<=basis.rangeJ(i).lastIndex(); ++k) {
                lawa::XType type = lawa::XWavelet;
                lawa::Index1D lambda(i, k, type);
                ++count;
                u(lambda, l) = count;
            }
        }
        printf("Adding wavelet coefficients\n");
        ht.addCoeff(u);
        std::cout << ht.getFrame() << std::endl;

        printf("\n...Printing Active Coefficients...\n");
        ht.printActive();
    }
    testIterators(ht);
}

void
apply(HTCoeffNode& u, const T& eps)
{
    T c = 0;
    HeHoOp op(basis, c);

    int numCols = u.getRank();
    HTCoeffFrame U(numCols), ret(numCols);

    typedef typename HTCoeffNode::const_iterator const_it;

    for (int j=1; j<=numCols; ++j) {
        for (const_it it=u.cbegin(j); it!=u.cend(j); ++it) {
            U((*it), j) = u((*it), j);
        }
        op.apply(U[j], eps, ret[j]);
    }

    u.setCoeff(ret);
    /*printf("Result frame is\n");
    for (int l=1; l<=numCols; ++l) {
        printf("Column %d\n", l);
        std::cout << ret[l] << std::endl;
    }*/
}


void
addCoefficients(Coefficients& w_p, const int& p,
                lawa::Coefficients<lawa::Bucket,T,lawa::Index1D>& v_bucket,
                const HTCoeffNode& U, const int& col)
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
apply_i(HTCoeffNode& U, const T& eps, const Coefficients& v)
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
    int numCols = U.getRank();
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
buildCoeff(Coefficients& v, const HTCoeffNode& U, const int& col)
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

    int count(0);
    T   scale(30);
    int j0 = basis.j0;
    bool flag(false);
    int offset(0);
    for (int l=1; l<=rank; ++l) {
        //if (l==3) scale = 0.5;
        if (l>3) scale = 40;
        for (int k=basis.mra.rangeI(j0).firstIndex();
                 k<=basis.mra.rangeI(j0).lastIndex(); ++k) {
            offset =    basis.mra.rangeI(j0).lastIndex()-2*k
                        +basis.mra.rangeI(j0).firstIndex();
            lawa::XType type = lawa::XBSpline;
            if (l==1) offset = 0;
            lawa::Index1D lambda(j0, k+offset, type);
            ++count;
            double r = ((double) rand() / (double) RAND_MAX);
            u(lambda, l) = r;//count/scale;
        }
        for (int i=basis.j0; i<=j; ++i) {
            for (int k=basis.rangeJ(i).firstIndex();
                     k<=basis.rangeJ(i).lastIndex(); ++k) {
                offset =    basis.rangeJ(i).lastIndex()-2*k
                            +basis.rangeJ(i).firstIndex();
                lawa::XType type = lawa::XWavelet;
                if (l==1) offset = 0;
                lawa::Index1D lambda(i, k+offset, type);
                ++count;
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
    for (int i=1; i<=copy.getRank(); ++i) {
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
    //printBasis(basis, 0);
    //testHTCoeffNode(0, 3);
    testApply(0, 3);
    return 0;
}
