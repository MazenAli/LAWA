#include <lawa/lawa.h>


typedef     double                              T;
typedef     lawa::Index2D                       Index;
typedef     lawa::IndexSet<Index>               IndexSet;
typedef     lawa::Coefficients<lawa::Lexicographical,
            T, Index>                           Coefficients;
typedef     lawa::Basis<T, lawa::Orthogonal,
            lawa::Interval, lawa::Multi>        Basis_XY;
//typedef     lawa::Basis<T, lawa::Primal,
//            lawa::Interval, lawa::Dijkema>      Basis_XY;
typedef     lawa::TensorBasis2D<lawa::Adaptive,
            Basis_XY, Basis_XY>
                                                Basis2D;
typedef     lawa::SeparableRHS2D<T, Basis2D>    RHSINTEGRAL;
typedef     lawa::SmoothRHSWithAlignedSing2D
            <T, Basis2D, lawa::SparseGridGP>    RHSINTEGRAL2;
typedef     lawa::NoPreconditioner<T, Index>    Preconditioner;
typedef     flens::DenseVector<flens::Array<T>>
                                                DenseVectorT;
typedef     flens::GeMatrix<flens::FullStorage
            <T, cxxblas::ColMajor>>             FullColMatrixT;
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
typedef     lawa::LocalOperator1D<Basis_XY,
            Basis_XY, RefLaplace1D, Laplace1D>  LOp_Lapl1D;
typedef     lawa::LocalOperator1D<Basis_XY,
            Basis_XY, RefIdentity1D,
            Identity1D>                         LOp_Id1D;
typedef     lawa::LocalOperator2D<LOp_Lapl1D,
            LOp_Id1D>                           LOp_Lapl_Id_2D;
typedef     lawa::LocalOperator2D<LOp_Id1D,
            LOp_Lapl1D>                         LOp_Id_Lapl_2D;
typedef     lawa::FlexibleCompoundLocalOperator
            <Index,
            lawa::AbstractLocalOperator2D<T>>   Flex_COp_2D;
typedef     lawa::H1NormPreconditioner2D<T,
            Basis2D>                            Prec2D;


double  a   = 1./35.;
double  b   = 1./35.;
double  mx  = 0.5;
double  my  = 0.5;


struct precon
{
    T
    operator() (const Index&)
    {
        return 1.;
    }
};


T
zero(T, T)
{
    return 0.;
}


double
ux(double x)
{
    //return std::exp(-1*std::pow(x-mx, 2.)/(a*a));
    return -1.*(x-0.5)*(x-0.5)+0.25;
    /*
    if (x >= 0.4 && x <= 0.6) {
        return 1.;
    } else {
        return 0.;
    }*/
}


double
uy(double y)
{
    //return std::exp(-1*std::pow(y-my, 2.)/(b*b));
    return -1.*(y-0.5)*(y-0.5)+0.25;
    /*
    if (y >= 0.4 && y <= 0.6) {
        return 1.;
    } else {
        return 0.;
    }*/
}


double
uxx(double x)
{
    return ux(x)*(std::pow(x-mx,2.)*4./(a*a*a*a)
           -2./(a*a));
}


double
uyy(double y)
{
    return uy(y)*(std::pow(y-my,2.)*4./(b*b*b*b)
           -2./(b*b));
}


double
fxy(double x, double y)
{
    //return -1.*(uxx(x)*uy(y)+uyy(y)*ux(x));
    return -2.*(std::pow(x-0.5,2.)+
              std::pow(y-0.5,2.)-
              0.5);
}


// Basis
int                                               d(4);
Basis_XY                                          basis_xy(d, 2);
Basis2D                                           basis2d(basis_xy, basis_xy);


// Input u
DenseVectorT                                      sing_support;
FullColMatrixT                                    nodeltas;
lawa::SeparableFunction2D<T>                      uxy(ux, sing_support,
                                                      uy, sing_support);
RHSINTEGRAL                                       _u(basis2d, uxy,
                                                       nodeltas, nodeltas, 100);
Preconditioner                                    prec;
lawa::RHS<T, Index, RHSINTEGRAL, Preconditioner>  u(_u, prec);


// Data f
lawa::Function2D<T>                                 __f(fxy, sing_support,
                                                          sing_support);
RHSINTEGRAL2                                        _f(basis2d, __f, 100);
Prec2D                                              P(basis2d);
lawa::RHS<T, Index, RHSINTEGRAL2, Preconditioner>   f(_f, prec);


// Operator A
Laplace1D                                         LaplaceBil(basis_xy);
RefLaplace1D                                      RefLaplaceBil(
                                                  basis_xy.refinementbasis);
Identity1D                                        IdentityBil(basis_xy);
RefIdentity1D                                     RefIdentityBil(
                                                  basis_xy.refinementbasis);
LOp_Lapl1D                                        lOp_Lapl1D(basis_xy,
                                                  basis_xy, RefLaplaceBil,
                                                  LaplaceBil);
LOp_Id1D                                          lOp_Id1D(basis_xy,
                                                  basis_xy, RefIdentityBil,
                                                  IdentityBil);
LOp_Lapl_Id_2D                                    localLaplId(lOp_Lapl1D,
                                                  lOp_Id1D);
LOp_Id_Lapl_2D                                    localIdLapl(lOp_Id1D,
                                                  lOp_Lapl1D);
std::vector<lawa::AbstractLocalOperator2D<T>*>    ops;


template <typename _Int, typename _Prec>
void
sample(IndexSet Lambda, lawa::RHS<T, Index, _Int, _Prec>& f,
       Coefficients& ret, T tol, T alpha = 0.7,
       std::size_t max_it = 1e+02)
{
    assert(Lambda.size() > 0);

    Coefficients    sweep;  // current extension sweep direction
    Coefficients    total;  // total extention residual

    // Initial step
    lawa::FillWithZeros(Lambda, sweep);
    lawa::FillWithZeros(Lambda, total);
    ret.clear();
    ret = f(Lambda);

    // Sampling loop
    for (std::size_t k = 0; k < max_it; ++k) {
        lawa::extendMultiTree(basis2d, sweep, total, "standard", true/*false*/, true);
        IndexSet diff = lawa::supp(total);

        for (auto& lambda : Lambda) {
            diff.erase(lambda);
        }

        Coefficients res = f(diff);
        T res_norm = res.norm((T)2.);

        if (res_norm <= tol) {
            std::cout << "Tolerance " << res_norm << " reached after "
                      << k+1 << " iterations\n";
            return;
        }

        // Bucket sort
        T P_Lambda = (T)0.;
        T thresh   = std::sqrt((T)1.-alpha*alpha)*res_norm /
                     std::sqrt(T(res.size()));
        lawa::Coefficients<lawa::Bucket, T, lawa::Index2D> buckets;
        buckets.bucketsort(res, thresh);

        Coefficients new_ind;
        for (std::size_t i = 0; i < buckets.bucket_ell2norms.size(); ++i) {
            P_Lambda += std::pow(buckets.bucket_ell2norms[i], 2.0L);
            buckets.addBucketToCoefficients(new_ind, i);

            if (P_Lambda >= std::pow(alpha*res_norm, (T)2.)) break;
        }

        // New index set
        IndexSet new_sweep;

        for (auto& coeff : new_ind) {
            if (ret.find(coeff.first) == ret.end()) {
                completeMultiTree(basis2d, coeff.first, ret,
                                  new_sweep, 0, true);
            }
        }

        #ifdef DEBUG_SAMPLER
            // Debug***
            std::cout << "\n\n*****Debug f sampler*****\n\n";
            std::cout << "Current iteration k = " << k+1 << std::endl;
            std::cout << "Current res_norm = " << res_norm << std::endl;
            std::cout << "Current ret_norm = " << ret.norm(2.) << std::endl;
            std::cout << "Current Lambda =>\n";
            std::cout << Lambda;
            std::cout << "Size Lambda = " << Lambda.size() << std::endl;
            std::cout << "Current sweep =>\n";
            std::cout << sweep;
            std::cout << "Size sweep = " << sweep.size() << std::endl;
            std::cout << "Current total =>\n";
            std::cout << total;
            std::cout << "Size total = " << total.size() << std::endl;
            std::cout << "Current diff =>\n";
            std::cout << diff;
            std::cout << "Size diff = " << diff.size() << std::endl;
            std::cout << "Current ret =>\n";
            std::cout << ret;
            std::cout << "Size ret = " << ret.size() << std::endl;
            std::cout << "End of iteration\n\n";
            // End debug***
        #endif

        sweep.clear();
        lawa::FillWithZeros(new_sweep, sweep);
        ret += f(new_sweep);
        Lambda = lawa::supp(ret);
    }

    std::cerr << "Warning! Max iterations " << max_it << " reached!\n";
}


void
eval(IndexSet Lambda, Flex_COp_2D& A,
     Coefficients& u, Coefficients& ret,
     T tol, T alpha = 0.7,
     std::size_t max_it = 1e+02)
{
    assert(Lambda.size() > 0);

    Coefficients    sweep;  // current extension sweep direction
    Coefficients    total;  // total extention residual

    // Initial step
    lawa::FillWithZeros(Lambda, sweep);
    lawa::FillWithZeros(Lambda, total);
    ret.clear();
    lawa::FillWithZeros(Lambda, ret);
    A.eval(u, ret, P);

    // Sampling loop
    for (std::size_t k = 0; k < max_it; ++k) {
        lawa::extendMultiTree(basis2d, sweep, total, "standard", true/*false*/, true);
        Coefficients res = total;
        A.eval(u, res, P);

        for (auto& lambda : Lambda) {
            res.erase(lambda);
        }

        T res_norm = res.norm((T)2.);

        if (res_norm <= tol) {
            std::cout << "Tolerance " << res_norm << " reached after "
                      << k+1 << " iterations\n";
            return;
        }

        // Bucket sort
        T P_Lambda = (T)0.;
        T thresh   = std::sqrt((T)1.-alpha*alpha)*res_norm /
                     std::sqrt(T(res.size()));
        lawa::Coefficients<lawa::Bucket, T, lawa::Index2D> buckets;
        buckets.bucketsort(res, thresh);

        Coefficients new_ind;
        for (std::size_t i = 0; i < buckets.bucket_ell2norms.size(); ++i) {
            P_Lambda += std::pow(buckets.bucket_ell2norms[i], 2.0L);
            buckets.addBucketToCoefficients(new_ind, i);

            if (P_Lambda >= std::pow(alpha*res_norm, (T)2.)) break;
        }

        // New index set
        IndexSet new_sweep;

        for (auto& coeff : new_ind) {
            if (ret.find(coeff.first) == ret.end()) {
                completeMultiTree(basis2d, coeff.first, ret,
                                  new_sweep, 0, true);
            }
        }

        #ifdef DEBUG_SAMPLER
            // Debug***
            std::cout << "\n\n*****Debug Au sampler*****\n\n";
            std::cout << "Current iteration k = " << k+1 << std::endl;
            std::cout << "Current res_norm = " << res_norm << std::endl;
            std::cout << "Current ret_norm = " << ret.norm(2.) << std::endl;
            std::cout << "Current Lambda =>\n";
            std::cout << Lambda;
            std::cout << "Size Lambda = " << Lambda.size() << std::endl;
            std::cout << "Current sweep =>\n";
            std::cout << sweep;
            std::cout << "Size sweep = " << sweep.size() << std::endl;
            std::cout << "Current total =>\n";
            std::cout << total;
            std::cout << "Size total = " << total.size() << std::endl;
            std::cout << "Current res =>\n";
            std::cout << res;
            std::cout << "Size res = " << res.size() << std::endl;
            std::cout << "Current ret =>\n";
            std::cout << ret;
            std::cout << "Size ret = " << ret.size() << std::endl;
            std::cout << "End of iteration\n\n";
            // End debug***
        #endif

        sweep.clear();
        lawa::FillWithZeros(new_sweep, sweep);
        ret.setToZero();
        A.eval(u, ret, P);
        Lambda = lawa::supp(ret);
    } // sampling loop

    std::cerr << "Warning! Max iterations " << max_it << " reached!\n";
}


int
main()
{
    basis_xy.enforceBoundaryCondition<lawa::DirichletBC>();
    Coefficients F, U, AU, res;
    IndexSet Lambda;
    T gamma     = 0.2;
    T tol_U     = 1e-08;
    T tol_F     = 1e-08;
    T tol_AU    = 1e-08;
    lawa::getSparseGridIndexSet(basis2d, Lambda, 1, 0, gamma);

    ops.push_back(&localLaplId);
    ops.push_back(&localIdLapl);
    Flex_COp_2D                                       A(ops);

    lawa::sample_f(basis2d, Lambda, u, U, tol_U, true);
    precon _p;
    lawa::saveCoeffVector2D(U, basis2d, "u.dat");
    lawa::plot2D(basis2d, U, _p, zero, 0., 1., 0., 1., 1e-02, "plot_u");
    lawa::sample_Au(basis2d, Lambda, A, P, U, AU, tol_AU, true);
    lawa::sample_f(basis2d, Lambda, f, P, F, tol_F, true);
    res = F - AU;
    for (auto& lambda : AU) {
        lambda.second /= P(lambda.first);
    }
    for (auto& lambda : F) {
        lambda.second /= P(lambda.first);
    }
    std::cout << "And the residual iiiiis... "
              << res.norm(2.) << std::endl;

    lawa::saveCoeffVector2D(AU, basis2d, "Au.dat");
    lawa::saveCoeffVector2D(F, basis2d, "f.dat");
    lawa::plot2D(basis2d, AU, _p, zero, 0., 1., 0., 1., 1e-02, "plot_Au");
    lawa::plot2D(basis2d, F, _p, zero, 0., 1., 0., 1., 1e-02, "plot_f");
    std::cout << "Everything exited fine\n";
    return 0;
}
