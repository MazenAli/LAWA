#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_HTAWGM_TCC
#define LAWA_METHODS_ADAPTIVE_SOLVERS_HTAWGM_TCC 1

#include <cassert>
#include <iostream>
#include <cmath>

#include <flens/flens.cxx>

#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/datastructures/sepcoefficients.h>
#include <lawa/methods/adaptive/algorithms/coeffops.h>
#include <lawa/methods/adaptive/algorithms/indexops.h>
#include <lawa/methods/adaptive/algorithms/multitreeoperations.h>

#ifndef MAX
    #define MAX(x, y) (x>y) ? x : y
#endif

namespace lawa
{

template <typename Optype, typename Basis, typename T>
unsigned
galerkin_pcg(Sepop<Optype>& A,
             Sepdiagscal<Basis>& S,
             HTCoefficients<T, Basis>& x,
             const HTCoefficients<T, Basis>& b,
             const std::vector<IndexSet<Index1D> >& Lambda,
             T& residual,
             const bool uzero,
             const T tol,
             const unsigned maxit,
             const T delta1,
             const T delta2,
             const T delta3,
             const T trunc)
{
    assert(A.dim()==S.dim());
    assert(A.dim()==(unsigned) x.dim());
    assert(A.dim()==(unsigned) b.dim());
    assert(A.dim()==Lambda.size());

    HTCoefficients<T, Basis>             r(b);
    HTCoefficients<T, Basis>             p(x.dim(),   x.basis(), x.map());
    p.tree().set_tree(x.tree());
    HTCoefficients<T, Basis>             tmp(x.dim(), x.basis(), x.map());
    tmp.tree().set_tree(x.tree());
    HTCoefficients<T, Basis>             rold(x.dim(), x.basis(), x.map());
    rold.tree().set_tree(x.tree());

    T nrmp, trunc_acc, nrmz, trunc_prec_s, t, nrmr, zr, trunc_prec;

    /* Initial residual */
    if (!uzero) {
        tmp = eval(A, x, Lambda, Lambda);
        scal(-1., tmp);
        r.tree() = add_truncate(tmp.tree(), r.tree(), trunc);
    }

    /* Single precondition truncation accuracy */
    t            = compOmegamax2(Lambda, S.order())/
                   compOmegamin2(S.basis(), S.dim(), S.order());
    trunc_prec_s = 10.*S.eps()/std::sqrt(t);
    #ifdef VERBOSE
        std::cout << "galerkin_pcg: trunc_prec single = " << trunc_prec_s
                  << std::endl;
    #endif
    nrmr      = nrm2(r);
    rold      = r;
    trunc_prec= std::min(1e-03, trunc_prec_s*nrmr);
    r         = eval(S, r, Lambda, trunc_prec);
    residual  = nrm2(r);
    trunc_acc = residual;

    p         = eval(S, r, Lambda, delta1*residual);
    zr        = dot(rold, p);
    nrmp      = nrm2(p);
    nrmz      = nrmp;

    T bk           = 1.;
    T trunc_search = trunc_acc;
    for (unsigned k=1; k<=maxit; ++k) {
        T ak, pAp;

        #ifdef VERBOSE
            std::cout << "galerkin_pcg: Iteration " << k
                      << " residual " << residual << std::endl;
            std::cout << "galerkin_pcg: max rank z "
                      << r.tree().max_rank() << std::endl;
            std::cout << "galerkin_pcg: max rank solution "
                      << x.tree().max_rank() << std::endl;
            std::cout << "galerkin_pcg: max rank p "
                      << p.tree().max_rank() << std::endl;
            std::cout << "galerkin_pcg: max rank b "
                      << b.tree().max_rank() << std::endl;
        #endif

        if (residual<=tol && k>1) {
            #ifdef VERBOSE
                std::cout << "galerkin_pcg: Tolerance reached r = "
                          << residual << std::endl;
            #endif
            return k;
        }

        /* alpha_k */
        tmp        = eval(A, p, Lambda, Lambda);
        ak         = dot(b, p);
        pAp        = dot(p, tmp);
        if (!uzero || k>1) {
            ak -= dot(x, tmp);
        }
        ak       /= pAp;
        trunc_acc = nrmz*delta3*std::fabs(ak);

        /* Update x_k and r_k */
        if (uzero && k==1) {
            x.tree() = ak*p.tree();
        } else {
            tmp = p;
            scal(ak, tmp);
            x.tree() = add_truncate(tmp.tree(), x.tree(), trunc_acc);
        }

        #ifdef VERBOSE
            std::cout << "galerkin_pcg: trunc_acc    = " << trunc_acc
                      << std::endl;
            std::cout << "galerkin_pcg: trunc_search = " << trunc_search
                      << std::endl;
            std::cout << "galerkin_pcg: nrmp         = " << nrmp
                      << std::endl;
            std::cout << "galerkin_pcg: nrmz         = " << nrmz
                      << std::endl;
            std::cout << "galerkin_pcg: ak           = " << ak
                      << std::endl;
            std::cout << "galerkin_pcg: bk           = " << bk
                      << std::endl;
            std::cout << "galerkin_pcg: nrmr         = " << nrmr
                      << std::endl;
        #endif

        r = eval(A, x, Lambda, Lambda);
        r.truncate(trunc_acc);
        scal(-1., r);
        scal(-1., rold);
        r.tree()   = add_truncate(b.tree(), r.tree(), trunc_acc);
        tmp.tree() = add_truncate(r.tree(), rold.tree(), trunc_acc);
        rold       = r;

        nrmr       = nrm2(r);
        trunc_prec = std::min(trunc_acc, trunc_prec_s*nrmr);
        r          = eval(S, r, Lambda, trunc_prec);
        residual   = nrm2(r);
        trunc_prec = std::min(trunc_acc, trunc_prec_s*residual);
        r          = eval(S, r, Lambda, trunc_prec);

        /* Update p_k */
        nrmz = nrm2(r);
        trunc_search = delta2*(nrmz*nrmz)/(std::fabs(bk)*nrmp);
        if (k==1) {
            trunc_search = delta1*nrmz;
        } else {
            trunc_search = std::min(trunc_search, delta1*nrmz);
        }

        bk = dot(r, tmp)/zr;
        zr = dot(r, rold);
        if (bk<=0) {
            p = r;
        } else {
            scal(bk, p);
            p.tree()  = add_truncate(r.tree(), p.tree(), trunc_search);
        }
        nrmp = nrm2(p);
    }

    std::cerr << "galerkin_pcg: Max iterations reached: maxit " << maxit
              << " r = " << residual << std::endl;
    return maxit;
}


template <typename Optype, typename Basis, typename T>
std::vector<IndexSet<Index1D> >
presidual(Sepop<Optype>& A,
          Sepdiagscal<Basis>& S,
                HTCoefficients<T, Basis>& u,
                HTCoefficients<T, Basis>& f,
                SepCoefficients<Lexicographical, T, Index1D>& fcp,
                HTCoefficients<T, Basis>& r,
          const SeparableRHSD<T, Basis>& fint,
          const std::vector<IndexSet<Index1D> >& current,
          const std::vector<IndexSet<Index1D> >& sweep,
                std::vector<IndexSet<Index1D> >& total,
          const T trunc)
{
    assert(A.dim()==S.dim());
    assert(A.dim()==(unsigned) u.dim());
    assert(A.dim()==(unsigned) f.dim());
    assert(A.dim()==(unsigned) r.dim());
    assert(A.dim()==current.size());
    assert(A.dim()==sweep.size());
    assert(A.dim()==total.size());

    typedef typename std::vector<IndexSet<Index1D> >::size_type size_type;

    std::vector<IndexSet<Index1D> > diff(current.size());
    HTCoefficients<T, Basis>        Au(u.dim(),    u.basis(), u.map());
    Au.tree().set_tree(u.tree());
    HTCoefficients<T, Basis>        prec(u.dim(),  u.basis(), u.map());
    prec.tree().set_tree(u.tree());

    /* Determine extended index set for evaluation */
    for (size_type j=0; j<current.size(); ++j) {
        extendMultiTree(u.basis(), sweep[j], total[j],
                        "standard", false);
        diff[j] = total[j];

        for (auto& lambda : current[j]) {
            diff[j].erase(lambda);
        }
    }

    /* Evaluate */
    genAddCoefficients(fcp, fint, diff);
    set(f, fcp);

    /* Compute residual */
    r = eval(A, u, total, current);
    scal(-1., r);
    r.tree() = add_truncate(f.tree(), r.tree(), trunc);
    r        = eval(S, r, total, trunc);

    return diff;
}


template <typename T, typename Basis>
std::vector<IndexSet<Index1D> >
bulk(const T alpha, const T resex,
           HTCoefficients<T, Basis>& res,
           std::vector<IndexSet<Index1D> >& Lambda,
     const std::vector<IndexSet<Index1D> >& diff)
{
    assert(Lambda.size()==(unsigned) res.dim());

    typedef typename std::vector<IndexSet<Index1D> >::size_type size_type;

    std::vector<flens::DenseVector<flens::Array<T> > >
                                    sigmas(Lambda.size());
    std::vector<Coefficients<Lexicographical, T, Index1D> >
                                    cont(Lambda.size());
    std::vector<IndexSet<Index1D> > sweep(Lambda.size());

    T thresh = std::sqrt((T) 1-alpha*alpha)*resex;

    /* Compute contractions */
    res.orthogonalize_svd(sigmas);
    contraction(res, diff, sigmas, cont);

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
    size_type num = 0;
    buckets.bucketsort(contvec, thresh);

    for (int i=buckets.bucket_ell2norms.size()-1; i>=0;
         --i) {
        P_Lambda += std::pow(buckets.bucket_ell2norms[i], 2.0L);
        if (P_Lambda > thresh*thresh) {
            num = i;
            break;
        }
    }

    for (size_type i=0; i<=num; ++i) {
        buckets.addBucketToCoefficients(newind, i);
    }

    /* Update index set and set new sweep */
    for (auto& it : newind) {
        completeMultiTree(res.basis(), it.first,
                          Lambda[it.first.d-1],
                          sweep[it.first.d-1],
                          true);
    }

    #ifdef VERBOSE
        restrict(res, Lambda);
        std::cout << "bulk: true alpha = " << nrm2(res)/resex << std::endl;
    #endif

    return sweep;
}


template <typename Optype, typename Basis, typename T>
unsigned
htawgm(Sepop<Optype>&                   A,
       Sepdiagscal<Basis>&              S,
       HTCoefficients<T, Basis>&        u,
       const SeparableRHSD<T, Basis>&   f,
       std::vector<IndexSet<Index1D> >& Lambda,
       double&                          residual,
       HTAWGM_Params&                   params)
{
    assert(A.dim()==S.dim());
    assert(A.dim()==(unsigned) u.dim());
    assert(A.dim()==Lambda.size());

    typedef typename std::vector<IndexSet<Index1D> >::size_type size_type;

    std::vector<IndexSet<Index1D> >              sweep(A.dim());
    std::vector<IndexSet<Index1D> >              total(A.dim());
    SepCoefficients<Lexicographical, T, Index1D> Fcp(f.rank(), f.dim());
    HTCoefficients<T, Basis>                     F(u.dim(),    u.basis(),
                                                   u.map());
    F.tree().set_tree(u.tree());
    HTCoefficients<T, Basis>                     r(u.dim(),    u.basis(),
                                                   u.map());
    r.tree().set_tree(u.tree());
    HTCoefficients<T, Basis>                     Au(u.dim(),   u.basis(),
                                                   u.map());
    Au.tree().set_tree(u.tree());

    sweep = Lambda;
    total = Lambda;

    genCoefficients(Fcp, f, Lambda);
    set(F, Fcp);
    r.tree() = F.tree();

    T nrmf = nrm2(F);

    /* Reduction parameter */
    T kappap = std::sqrt(2.*(T)S.dim()-3.);
    T alpha  = params.recompr;
    T kappar = params.theta/(1.+kappap*(1.+alpha));

    #ifdef VERBOSE
        std::cout << "htawgm: kappar = " << kappar << std::endl;
    #endif

    /* Initial residual */
    if (!params.uzero) {
        Au = eval(A, u, Lambda, Lambda);
        Au.truncate(nrmf*1e-01);
        scal(-1., Au);
        r.tree() = add_truncate(F.tree(), Au.tree(), nrmf*1e-01);
    }

    r         = eval(S, r, Lambda, params.tol_awgm*1e-01);
    residual  = nrm2(r);

    if (residual<=params.tol_awgm) {
        #ifdef VERBOSE
            std::cout << "htawgm: Tolerance reached r = " << residual
                      << std::endl;
        #endif

        return 0;
    }

    T tol;
    for (unsigned k=1; k<=params.maxit_awgm; ++k) {
        unsigned pcg_it, size;
        T        res_pcg;

        #ifdef VERBOSE
            std::cout << "htawgm: Iteration " << k
                      << " r = " << residual << std::endl;
        #endif

        /* Galerkin solve */
        tol    = params.gamma*residual;
        pcg_it = galerkin_pcg(A, S, u, F, Lambda, res_pcg,
                              params.uzero,
                              tol,
                              params.maxit_pcg,
                              params.delta1_pcg,
                              params.delta2_pcg,
                              params.delta3_pcg,
                              tol*1e-01);
        #ifdef VERBOSE
            std::cout << "htawgm: galerkin_pcg required " << pcg_it
                      << " iterations to reach tolerance "
                      << tol
                      << std::endl;
        #else
            (void) pcg_it;
        #endif

        /* Approximate residual */
        sweep        = presidual(A, S, u, F, Fcp, r, f,
                                 Lambda, sweep, total,
                                 params.tol_awgm);
        residual     = nrm2(r);

        if (residual<=params.tol_awgm) {
            #ifdef VERBOSE
                std::cout << "htawgm: Tolerance reached r = " << residual
                          << std::endl;
            #endif

            return k;
        }

        /* Bulk chasing */
        sweep = bulk(params.alpha, residual,
                     r, Lambda, sweep);

        /* New RHS */
        restrict(F, Lambda);

        /* Extend u to new Lambda */
        extend(u, Lambda);
        params.uzero = false;

        #ifdef VERBOSE
            std::cout << "htawgm: Index set sizes\n";
            size = 0;
            for (size_type j=0; j<Lambda.size(); ++j) {
                std::cout << "htawgm: d = " << j+1
                          << " : " << Lambda[j].size() << std::endl;
                size += Lambda[j].size();
                FLENS_DEFAULT_INDEXTYPE jmax = 0;
                for (const auto& it : Lambda[j]) {
                    FLENS_DEFAULT_INDEXTYPE level = it.j;
                    if (it.xtype==XWavelet) ++level;
                    jmax = MAX(level, jmax);
                }
                std::cout << "htawgm: jmax = " << jmax << std::endl;
            }
            std::cout << "htawgm: Overall = " << size << std::endl;
        #else
            (void) size;
        #endif
    }

    std::cerr << "htawgm: Max iterations reached: maxit "
              << params.maxit_awgm << " r = " << residual << std::endl;

    return params.maxit_awgm;
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_HTAWGM_TCC 1
