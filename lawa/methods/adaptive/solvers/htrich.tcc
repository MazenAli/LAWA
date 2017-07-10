#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_HTRICH_TCC
#define LAWA_METHODS_ADAPTIVE_SOLVERS_HTRICH_TCC 1

#include <flens/flens.cxx>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/algorithms/coeffops.h>
#include <lawa/methods/adaptive/algorithms/indexops.h>
#include <lawa/methods/adaptive/algorithms/multitreeoperations.h>

#include <chrono>


namespace lawa
{

/* Implementation according to Dahmen/Bachmayr */
template <typename Optype, typename Basis, typename T>
unsigned
htrich(      Sepop<Optype>&                   A,
             Sepdiagscal<Basis>&              S,
             HTCoefficients<T, Basis>&        u,
             SeparableRHSD<T, Basis>&         f,
             std::vector<IndexSet<Index1D> >& Lambda,
             T&                               residual,
       const HTRICH_Params&                   params)
{
    assert(A.dim()==S.dim());
    assert(A.dim()==(unsigned) u.dim());
    assert(A.dim()==Lambda.size());

    typedef typename std::vector<IndexSet<Index1D> >::size_type size_type;

    SepCoefficients<Lexicographical, T, Index1D> Fcp(f.rank(), f.dim());
    HTCoefficients<T, Basis>                     F(u.dim(),    u.basis(),
                                                   u.map());
    F.tree().set_tree(u.tree());
    HTCoefficients<T, Basis>                     r(u.dim(),    u.basis(),
                                                   u.map());
    r.tree().set_tree(u.tree());
    HTCoefficients<T, Basis>                     Au(u.dim(),   u.basis(),
                                                   u.map());
    std::vector<IndexSet<Index1D> >              total(Lambda.size());
    Au.tree().set_tree(u.tree());

    genCoefficients(Fcp, f, Lambda);
    set(F, Fcp);

    T eta;
    T cv     = (1.-S.eps())/(4.*params.nrmA);
    auto SF  = applyScale(S, F, Lambda, 1e-04);
    T nrmf   = nrm2(SF);
    auto J   = findminj(params.rho, params.omega, params.beta1,
                        params.beta2, params.kappa1, params.maxit_inner);
    auto start = std::chrono::system_clock::now();
    for (unsigned k=0; k<params.maxit_rich; ++k) {
        T thresh = params.kappa1*params.eps0*std::pow(2., -1.*k-1.);

        /* Inner Richardson iterations */
        for (unsigned j=0; j<J; ++j) {
            T etakj = std::pow(params.rho, j+1)*
                      std::pow(2., -1.*k)*params.eps0;

            /* Evaluate residual */
            if (k==0 && j==0 && params.uzero) {
                eta      = .5*etakj*cv/nrmf;
                S.set_nu(eta);
                T tol    = .5*eta;
                r        = applyScale(S, F, Lambda, tol);
                residual = nrm2(r);
            } else {
                T nrmu   = nrm2(u);
                T max    = std::max(nrmu, nrmf);
                eta      = .5*etakj*cv/max;
                S.set_nu(eta);
                T tol    = .5*eta;
                presidual(A, S, u, F, Fcp, r, f, Lambda,
                           total, tol);
                residual = nrm2(r);
            }

            auto end     = std::chrono::system_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start);
            std::cout << "htrich: Outer iteration " << k
                      << ", inner iteration " << j
                      << ", current elapsed time "
                      << elapsed.count() << " seconds\n";

            #ifdef VERBOSE
                std::cout << "htrich: Outer iteration " << k
                          << ", inner iteration " << j
                          << ", residual "
                          << residual << std::endl;
            #endif

            /* Update solution */
            std::cout << "etakj => " << etakj << std::endl;
            scal(params.omega, r);
            if (k==0 && j==0 && params.uzero) {
                u = r;
            } else {
                Lambda    = total;
                extend(u, Lambda);
                T tol     = params.beta1*etakj;
                u.tree()  = add_truncate(u.tree(), r.tree(), tol);
                tol       = params.beta2*etakj;
                Lambda    = coarsen(u, Lambda, tol);
                restrict(u, Lambda);
                total     = Lambda;
            }

            /* Display new index sets and ranks */
            #ifdef VERBOSE
                std::cout << "htrich: max rank solution "
                          << u.tree().max_rank() << std::endl;
                std::cout << "htrich: Index set sizes\n";
                unsigned size = 0;
                for (size_type j=0; j<Lambda.size(); ++j) {
                    std::cout << "htrich: d = " << j+1
                              << " : " << Lambda[j].size() << std::endl;
                    size += Lambda[j].size();
                    FLENS_DEFAULT_INDEXTYPE jmax = 0;
                    for (const auto& it : Lambda[j]) {
                        FLENS_DEFAULT_INDEXTYPE level = it.j;
                        if (it.xtype==XWavelet) ++level;
                        jmax = MAX(level, jmax);
                    }
                    std::cout << "htrich: jmax = " << jmax << std::endl;
                }
                std::cout << "htrich: Overall = " << size << std::endl;
            #endif

            T test = params.cA*params.rho*residual+
                     etakj*(params.cA*params.rho+params.omega+
                            params.beta1+params.beta2);
            if (test<=thresh) {
                break;
            }
        }

        /* Here follows coarse/compress */
        T tol     = params.kappa2*params.eps0*std::pow(2., -1.*k-1.);
        std::cout << "Recompress tol " << tol << std::endl;
        u.truncate(tol);
        tol       = params.kappa3*params.eps0*std::pow(2., -1.*k-1.);
        std::cout << "Coarse tol " << tol << std::endl;
        Lambda    = coarsen(u, Lambda, tol);
        restrict(u, Lambda);
        total     = Lambda;

        /* Check error */
        T test = std::pow(2., -k)*params.eps0;
        if (test<=params.tol_rich) return k;
    }

    std::cerr << "htrich: Max iterations reached: maxit "
              << params.maxit_rich << " r = " << residual << std::endl;

    return params.maxit_rich;
}

template <typename T, typename Basis>
std::vector<IndexSet<Index1D> >
coarsen(      HTCoefficients<T, Basis>&         u,
        const std::vector<IndexSet<Index1D> >&  Lambda,
        const T                                 eps)
{
    assert(Lambda.size()==(unsigned) u.dim());

    typedef typename std::vector<IndexSet<Index1D> >::size_type size_type;

    std::vector<flens::DenseVector<flens::Array<T> > >
                                    sigmas(Lambda.size());
    std::vector<Coefficients<Lexicographical, T, Index1D> >
                                    cont(Lambda.size());
    std::vector<IndexSet<Index1D> > ret(Lambda.size());

    T thresh = eps*eps;

    /* Compute contractions */
    u.orthogonalize_svd(sigmas);
    contraction(u, Lambda, sigmas, cont);

    /* Distribute contractions */
    Coefficients<Lexicographical, T, Index1DC> contvec;
    for (size_type j=0; j<cont.size(); ++j) {
        for (auto& lambda : cont[j]) {
            Index1DC index(lambda.first.j,
                           lambda.first.k,
                           lambda.first.xtype,
                           j+1);
            contvec[index] = lambda.second;
        }
    }

    /* Bucket sort */
    Coefficients<Bucket, T, Index1DC>           buckets;
    Coefficients<Lexicographical, T, Index1DC>  newind;
    T P_Lambda    = 0;
    int num       = 0;
    buckets.bucketsort(contvec, eps);

    for (int i=buckets.bucket_ell2norms.size()-1; i>=0;
         --i) {
        P_Lambda += std::pow(buckets.bucket_ell2norms[i], 2.0L);
        if (P_Lambda > thresh) {
            num = i;
            break;
        }
    }

    for (int i=0; i<=num; ++i) {
        buckets.addBucketToCoefficients(newind, i);
    }

//    /* Remove indices from last bucket */
//    Coefficients<Lexicographical, T, Index1DC>  remove;
//    buckets.addBucketToCoefficients(newind, num);
//    buckets.addBucketToCoefficients(remove, num);
//    for (auto& lambda : remove) {
//        P_Lambda += lambda.second*lambda.second;
//        if (P_Lambda>thresh) {
//            break;
//        }
//        if (count>1) newind.erase(lambda.first);
//        ++count;
//    }

    /* Set new index set */
    for (auto& it : newind) {
        completeMultiTree(u.basis(), it.first,
                          ret[it.first.d-1],
                          ret[it.first.d-1],
                          true);
    }

    return ret;
}

template <typename Optype, typename Basis, typename T>
void
presidual(Sepop<Optype>& A,
          Sepdiagscal<Basis>& S,
                HTCoefficients<T, Basis>& u,
                HTCoefficients<T, Basis>& f,
                SepCoefficients<Lexicographical, T, Index1D>& fcp,
                HTCoefficients<T, Basis>& r,
                SeparableRHSD<T, Basis>& fint,
          const std::vector<IndexSet<Index1D> >& current,
                std::vector<IndexSet<Index1D> >& total,
          const T trunc)
{
    assert(A.dim()==S.dim());
    assert(A.dim()==(unsigned) u.dim());
    assert(A.dim()==(unsigned) f.dim());
    assert(A.dim()==(unsigned) r.dim());
    assert(A.dim()==current.size());
    assert(A.dim()==total.size());

    typedef typename std::vector<IndexSet<Index1D> >::size_type size_type;

    std::vector<IndexSet<Index1D> > eval_diff(current.size());

    /* Determine extended index set for evaluation */
    for (size_type j=0; j<current.size(); ++j) {
        extendMultiTree(u.basis(), current[j], total[j],
                        "standard", false);
        eval_diff[j] = total[j];
        for (auto& lambda : fcp(1, j+1)) {
            eval_diff[j].erase(lambda.first);
        }
    }

    /* Evaluate */
    genAddCoefficients(fcp, fint, eval_diff);
    set(f, fcp, total);

    /* Compute residual */
    r = evaleff2(A, S, u, total, current, trunc);
    scal(-1., r);
    auto tmp = applyScale(S, f, total, trunc);
    r.tree() = add_truncate(tmp.tree(), r.tree(), trunc);
}

template <typename T>
unsigned
findminj(const T rho, const T omega, const T beta1,
         const T beta2, const T kappa1,
         unsigned max)
{
    unsigned j=0;

    for (; j<=max; ++j) {
        T test = std::pow(rho, (T) j)*(1.+(omega+beta1+beta2)*(T) j);
        if (test<=0.5*kappa1) return j;
    }

    return j;
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_HTRICH_TCC
