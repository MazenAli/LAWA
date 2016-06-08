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
             const T delta,
             const T trunc)
{
    assert(A.dim()==S.dim());
    assert(A.dim()==(unsigned) x.dim());
    assert(A.dim()==(unsigned) b.dim());
    assert(A.dim()==Lambda.size());

    HTCoefficients<T, Basis>             r(x.dim(), x.basis());
    HTCoefficients<T, Basis>             Sb(x.dim(), x.basis());
    HTCoefficients<T, Basis>             p(x.dim(), x.basis());
    HTCoefficients<T, Basis>             Ap(x.dim(), x.basis());
    HTCoefficients<T, Basis>             tmp(x.dim(), x.basis());

    T   nrmp;

    Sb = eval(S, b, Lambda, trunc);
    r  = Sb;

    /* Initial residual */
    if (!uzero) {
        Ap = eval(S, x, Lambda, trunc);
        Ap = eval(A, Ap, Lambda, Lambda);
        Ap.truncate(trunc);
        Ap = eval(S, Ap, Lambda, trunc);
        scal(-1., Ap);
        r.tree() = add_truncate(Ap.tree(), r.tree(), trunc);
    }

    residual = nrm2(r);
    p        = r;
    nrmp     = residual;

    for (unsigned k=1; k<=maxit; ++k) {
        T ak, bk, pAp, trunc_acc;

        #ifdef VERBOSE
            std::cout << "galerkin_pcg: Iteration " << k
                      << " residual " << residual << std::endl;
            std::cout << "galerkin_pcg: max rank r "
                      << r.tree().max_rank() << std::endl;
            std::cout << "galerkin_pcg: max rank solution "
                      << x.tree().max_rank() << std::endl;
            std::cout << "galerkin_pcg: max rank p "
                      << p.tree().max_rank() << std::endl;
            std::cout << "galerkin_pcg: max rank Sb "
                      << Sb.tree().max_rank() << std::endl;
        #endif

        if (residual<=tol) {
            #ifdef VERBOSE
                std::cout << "galerkin_pcg: Tolerance reached r = "
                          << residual << std::endl;
            #endif

            return k;
        }

        /* alpha_k */
        trunc_acc = delta*nrmp;
        Ap        = eval(S, p, Lambda, trunc_acc);
        std::cout << "galerkin_pcg: max rank Sp " << Ap.tree().max_rank() << std::endl;
        Ap        = eval(A, Ap, Lambda, Lambda);
        std::cout << "galerkin_pcg: max rank ASp " << Ap.tree().max_rank() << std::endl;
        Ap.truncate(MAX(trunc, trunc_acc));
        Ap        = eval(S, Ap, Lambda, trunc_acc);
        std::cout << "galerkin_pcg: max rank SASp " << Ap.tree().max_rank() << std::endl;
        ak        = dot(Sb, p);
        std::cout << "(Sb, p) = " << ak << std::endl;

        if (!uzero || k>1) {
            ak -= dot(x, Ap);
        }

        std::cout << "(Sb, p)-(x, Ap) = " << ak << std::endl;
        pAp        = dot(p, Ap);
        std::cout << "pAp = " << pAp << std::endl;
        ak        /= pAp;
        trunc_acc *= ak;

        std::cout << "truncation accuracy = " << trunc_acc << std::endl;
        std::cout << "nrmp  = " << nrmp << std::endl;
        std::cout << "delta = " << delta << std::endl;
        std::cout << "ak    = " << ak << std::endl;

        /* Update x_k and r_k */
        if (uzero && k==1) {
            x.tree() = ak*p.tree();
        } else {
            tmp = p;
            scal(ak, tmp);
            x.tree() = add_truncate(tmp.tree(), x.tree(), trunc_acc);
        std::cout << "galerkin_pcg: max rank new Sx " << x.tree().max_rank() << std::endl;
        }

        r = eval(S, x, Lambda, trunc_acc);
        std::cout << "galerkin_pcg: max rank new Sx " << r.tree().max_rank() << std::endl;
        r = eval(A, r, Lambda, Lambda);
        std::cout << "galerkin_pcg: max rank new ASx " << r.tree().max_rank() << std::endl;
        r.truncate(trunc_acc);
        r = eval(S, r, Lambda, trunc_acc);
        std::cout << "galerkin_pcg: max rank new SASx " << r.tree().max_rank() << std::endl;
        scal(-1., r);
        r.tree() = add_truncate(Sb.tree(), r.tree(), trunc_acc);
        std::cout << "galerkin_pcg: max rank new Sb-SASx " << r.tree().max_rank() << std::endl;
        residual = nrm2(r);

        /* Update p_k */
        bk        = -dot(r, Ap)/pAp;
        scal(bk, p);
        Ap.tree() = r.tree()+p.tree();
        std::cout << "galerkin_pcg: max rank p before trunc " << Ap.tree().max_rank() << std::endl;
        p.tree()  = add_truncate(r.tree(), p.tree(), delta*nrm2(Ap));
        nrmp      = nrm2(p);
    }

    std::cerr << "galerkin_pcg: Max iterations reached: maxit " << maxit
              << " r = " << residual << std::endl;

    return maxit;
}




template <typename Optype, typename Basis, typename T>
unsigned
galerkin_pcg2(Sepop<Optype>& A,
             Sepdiagscal<Basis>& S,
             HTCoefficients<T, Basis>& x,
             const HTCoefficients<T, Basis>& b,
             const std::vector<IndexSet<Index1D> >& Lambda,
             T& residual,
             const bool uzero,
             const T tol,
             const unsigned maxit,
             const T delta,
             const T trunc)
{
    assert(A.dim()==S.dim());
    assert(A.dim()==(unsigned) x.dim());
    assert(A.dim()==(unsigned) b.dim());
    assert(A.dim()==Lambda.size());

    HTCoefficients<T, Basis>             r(b);
    HTCoefficients<T, Basis>             p(x.dim(), x.basis());
    HTCoefficients<T, Basis>             Ap(x.dim(), x.basis());
    HTCoefficients<T, Basis>             tmp(x.dim(), x.basis());

    T   nrmp;

    T factor = compUnDistFac(r, S.order());

    /* Initial residual */
    if (!uzero) {
        Ap = eval(A, x, Lambda, Lambda);
        scal(-1., Ap);
        r.tree() = add_truncate(Ap.tree(), r.tree(), trunc);
    }

    residual  = nrm2(r);
    residual *= factor;
    p         = evalS2(S, r, Lambda, trunc);
    nrmp      = nrm2(p);


    for (unsigned k=1; k<=maxit; ++k) {
        T ak, bk, pAp, trunc_acc;

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

        if (residual<=tol) {
            #ifdef VERBOSE
                std::cout << "galerkin_pcg: Tolerance reached r = "
                          << residual << std::endl;
            #endif

            return k;
        }

        /* alpha_k */
        Ap         = eval(A, p, Lambda, Lambda);
        trunc_acc  = delta*nrmp;
        Ap.truncate(MAX(trunc, trunc_acc));
        ak         = dot(b, p);
        pAp        = dot(p, Ap);

        if (!uzero || k>1) {
            ak -= dot(x, Ap);
        }

        ak        /= pAp;
        trunc_acc *= ak;
        std::cout << "trunc_acc = " << trunc_acc << std::endl;
        std::cout << "nrmp      = " << nrmp <<std::endl;
        std::cout << "ak        = " << ak << std::endl;
        std::cout << "bk        = " << bk << std::endl;

        /* Update x_k and r_k */
        if (uzero && k==1) {
            x.tree() = ak*p.tree();
        } else {
            tmp = p;
            scal(ak, tmp);
            x.tree() = add_truncate(tmp.tree(), x.tree(), trunc_acc);
        }

        r = eval(A, x, Lambda, Lambda);
        r.truncate(trunc_acc);
        scal(-1., r);
        r.tree() = add_truncate(b.tree(), r.tree(), trunc_acc);
        residual = nrm2(r);
        residual *= factor;
        r = evalS2(S, r, Lambda, trunc_acc);

        /* Update p_k */
        bk        = -dot(r, Ap)/pAp;
        scal(bk, p);
        Ap.tree() = r.tree()+p.tree();
        p.tree()  = add_truncate(r.tree(), p.tree(), delta*nrm2(Ap));
        nrmp      = nrm2(p);
    }

    std::cerr << "galerkin_pcg: Max iterations reached: maxit " << maxit
              << " r = " << residual << std::endl;

    return maxit;
}


template <typename Optype, typename Basis, typename T>
std::vector<IndexSet<Index1D> >
presidual(Sepop<Optype>& A,
          Sepdiagscal<Basis>& S,
          const HTCoefficients<T, Basis>& u,
                HTCoefficients<T, Basis>& f,
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
    SepCoefficients<Lexicographical, T, Index1D>
                                    fdiff_(fint.rank(),
                                            fint.dim());
    HTCoefficients<T, Basis>        fdiff(u.dim(), u.basis());
    HTCoefficients<T, Basis>        Au(u.dim(), u.basis());
    HTCoefficients<T, Basis>        prec(u.dim(), u.basis());
    htucker::HTuckerTree<T>         gram;
    flens::DenseVector
          <flens::Array<T> >        eps2vec(u.tree().depth()-1);

    /* Determine extended index set for evaluation */
    for (size_type j=0; j<current.size(); ++j) {
        extendMultiTree(u.basis(), sweep[j], total[j],
                        "standard", true);
        diff[j] = total[j];

        for (auto& lambda : current[j]) {
            diff[j].erase(lambda);
        }
    }

    /* Evaluate on difference only */
    genCoefficients(fdiff_, fint, diff);
    set(fdiff, fdiff_);
    f.tree() = add_truncate(f.tree(), fdiff.tree(), trunc);

    /* Compute residual */
    r = eval(A, u, total, current);
    r.truncate(trunc);
    scal(-1., r);
    r.tree() = add_truncate(f.tree(), r.tree(), trunc);

    return diff;
}


template <typename T, typename Basis>
std::vector<IndexSet<Index1D> >
bulk(const T alpha, const T resex, const T rescg,
           HTCoefficients<T, Basis>& res,
           std::vector<IndexSet<Index1D> >& Lambda,
     const std::vector<IndexSet<Index1D> >& diff,
     const std::vector<IndexSet<Index1D> >& total,
     const unsigned maxit)
{
    assert(Lambda.size()==(unsigned) res.dim());

    typedef typename std::vector<IndexSet<Index1D> >::size_type size_type;

    HTCoefficients<T, Basis>        res_current(res);
    std::vector<flens::DenseVector<flens::Array<T> > >
                                    sigmas(Lambda.size());
    std::vector<Coefficients<Lexicographical, T, Index1D> >
                                    cont(Lambda.size());
    std::vector<IndexSet<Index1D> > sweep(Lambda.size());
    std::vector<IndexSet<Index1D> > diff_current(diff);
    T current_nrm, bulk, thresh;

    bulk        = alpha*resex;
    current_nrm = rescg;
    thresh      = std::sqrt((T) 1-alpha*alpha)*resex;

    for (unsigned k=1; k<=maxit; ++k) {

        /* Compute contractions */
        res_current.orthogonalize_svd(sigmas);
        contraction(res_current, diff_current, sigmas, cont);

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
        T P_Lambda;
        P_Lambda = current_nrm*current_nrm;
        buckets.bucketsort(contvec, thresh);

        for (size_type i=0; i<buckets.bucket_ell2norms.size(); ++i) {
            P_Lambda += std::pow(buckets.bucket_ell2norms[i], 2.0L);
            buckets.addBucketToCoefficients(newind, i);

            if (P_Lambda >= bulk*bulk) break;
        }

        /* Update index set and set new sweep */
        for (auto& it : newind) {
            completeMultiTree(res.basis(), it.first,
                              Lambda[it.first.d-1],
                              sweep[it.first.d-1],
                              true);
        }
        restrict(res_current, Lambda);
        current_nrm = nrm2(res_current);

        /* Check bulk criterion */
        if (current_nrm>=bulk) {
            #ifdef VERBOSE
                std::cout << "bulk: Required " << k
                          << " iterations" << std::endl;
            #endif
            restrict(res, Lambda);
            return sweep;
        }

        /* Prepare next iteration */
        res_current.tree() = res.tree();
        for (size_type j=0; j<Lambda.size(); ++j) {
            diff_current[j] = total[j];
            for (auto& lambda : Lambda[j]) {
                diff_current[j].erase(lambda);
            }
        }
    }

    std::cerr << "bulk: Max iterations reached: maxit " << maxit
              << ", alpha = "  << alpha
              << ", rnew/rext = " << current_nrm/resex << std::endl;
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
    HTCoefficients<T, Basis>                     F(u.dim(), u.basis());
    HTCoefficients<T, Basis>                     Fnew(u.dim(), u.basis());
    HTCoefficients<T, Basis>                     r(u.dim(), u.basis());
    HTCoefficients<T, Basis>                     Au(u.dim(), u.basis());
    HTCoefficients<T, Basis>                     prec(u.dim(), u.basis());
    T                                            factor;

    sweep = Lambda;
    total = Lambda;

    genCoefficients(Fcp, f, Lambda);
    set(F, Fcp);
    r.tree() = F.tree();
    r.truncate(params.trunc_pres);

    factor = compUnDistFac(F, S.order());

    /* Initial residual */
    if (!params.uzero) {
        Au = eval(A, u, Lambda, Lambda);
        Au.truncate(params.trunc_pres);
        scal(-1., Au);
        r.tree() = add_truncate(F.tree(), Au.tree(), params.trunc_pres);
    }

    residual  = nrm2(r);
    residual *= factor;

    if (residual<=params.tol_awgm) {
        #ifdef VERBOSE
            std::cout << "htawgm: Tolerance reached r = " << residual
                      << std::endl;
        #endif

        return 0;
    }

    for (unsigned k=1; k<=params.maxit_awgm; ++k) {
        unsigned pcg_it, size;
        T        res_pcg;

        #ifdef VERBOSE
            std::cout << "htawgm: Iteration " << k
                      << " r = " << residual << std::endl;
        #endif

        /* Galerkin solve */
        pcg_it = galerkin_pcg2(A, S, u, F, Lambda, res_pcg,
                              params.uzero,
                              params.gamma*residual,
                              params.maxit_pcg,
                              params.delta_pcg,
                              params.gamma*residual*1e-01);
        #ifdef VERBOSE
            std::cout << "galerkin_pcg required " << pcg_it
                      << " iterations to reach tolerance "
                      << params.gamma*residual
                      << std::endl;
        #else
            (void) pcg_it;
        #endif

        /* Approximate residual */
        sweep = presidual(A, S, u, F, r, f,
                          Lambda, sweep, total,
                          params.trunc_pres);

        residual  = nrm2(r);
        factor    = compUnDistFac(r, S.order());
        std::cout << "factor = " << factor << std::endl;
        residual *= factor;


        if (residual<=params.tol_awgm) {
            #ifdef VERBOSE
                std::cout << "htawgm: Tolerance reached r = " << residual
                          << std::endl;
            #endif

            return k;
        }

        /* Bulk chasing */
        sweep = bulk(params.alpha, residual, res_pcg,
                     r, Lambda, sweep, total,
                     params.maxit_bulk);

        /* New RHS */
        restrict(F, Lambda);
        F.truncate(params.trunc_pres);

        /* Extend u to new Lambda */
        extend(u, Lambda);
        params.uzero = false;

        #ifdef VERBOSE
            std::cout << "ht_awgm: Index set sizes\n";
            size = 0;
            for (size_type j=0; j<Lambda.size(); ++j) {
                std::cout << "d = " << j+1
                          << " : " << Lambda[j].size() << std::endl;
                size += Lambda[j].size();
            }
            std::cout << "Overall = " << size << std::endl;
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
