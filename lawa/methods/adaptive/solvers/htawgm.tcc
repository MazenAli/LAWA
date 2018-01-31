#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_HTAWGM_TCC
#define LAWA_METHODS_ADAPTIVE_SOLVERS_HTAWGM_TCC 1

#include <cassert>
#include <iostream>
#include <cmath>
#include <chrono>

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

    T nrmp, trunc_acc, nrmz, trunc_prec_s, t, nrmr, zr;

    /* Initial residual */
    if (!uzero) {
        tmp = eval(A, x, Lambda, Lambda);
        scal(-1., tmp);
        r.tree() = add_truncate(tmp.tree(), r.tree(), trunc);
    }

    t            = compOmegamax2(Lambda, S.order());
    trunc_prec_s = 1./std::sqrt(t);

    rold      = r;
    nrmr      = nrm2(r);
    r         = applyScale(S, r, Lambda, trunc*nrmr);
    residual  = nrm2(r);
    trunc_acc = trunc*nrmr;

    p         = applyScale(S, r, Lambda, trunc*residual);
    zr        = 1.;
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
        if (k==1) {
            trunc_acc = std::max(delta3*std::fabs(ak)*tol*1e-01,
                                 nrmz*delta3*std::fabs(ak));
        } else {
            trunc_acc = nrmz*delta3*std::fabs(ak);
        }

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
        T eps_     = nrmp*1e-01;
        r          = applyScale(S, r, Lambda, trunc_prec_s*eps_*0.5);
        residual   = nrm2(r);
        r          = applyScale(S, r, Lambda, trunc_search);

        /* Update p_k */
        nrmz = nrm2(r);
        trunc_search = delta2*(nrmz*nrmz)/(std::fabs(bk)*nrmp);
        if (k==1) {
            trunc_search = trunc_acc;
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
unsigned
galerkin_pcg2(      Sepop<Optype>& A,
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
    HTCoefficients<T, Basis>             p(x.dim(),   x.basis(), x.map());
    p.tree().set_tree(x.tree());
    HTCoefficients<T, Basis>             tmp(x.dim(), x.basis(), x.map());
    tmp.tree().set_tree(x.tree());

    T nrmp, trunc_acc;

    /* Initial residual */
    auto tgal0 = std::chrono::system_clock::now();
    auto stiff = assemble_laplace<T, Optype>(A, Lambda, Lambda);
    auto tgal1 = std::chrono::system_clock::now();
    auto dgal  = std::chrono::duration_cast<std::chrono::seconds>(tgal1-tgal0);
    std::cout << "pcg: assemble time: " << dgal.count() << std::endl;

    tgal0 = std::chrono::system_clock::now();

    if (!uzero) {
        tmp = evallaplace(stiff, S, x, Lambda, Lambda, trunc/2.);
        scal(-1., tmp);
        r.tree() = add_truncate(tmp.tree(), b.tree(), trunc/2.);
    }

    residual  = nrm2(r);
    tgal1 = std::chrono::system_clock::now();
    dgal  = std::chrono::duration_cast<std::chrono::seconds>(tgal1-tgal0);
    std::cout << "pcg: initial residual time: " << dgal.count() << std::endl;
    trunc_acc = residual*delta*1e-01;
    p         = r;
    nrmp      = residual;

    #ifdef VERBOSE
        std::cout << "pcg: Iteration " << 0
                  << " residual " << residual << std::endl;
        std::cout << "pcg: max rank r "
                  << r.tree().max_rank() << std::endl;
    #endif

    for (unsigned k=1; k<=maxit; ++k) {
        T ak, pAp, bk, rAp;
        std::vector<HTCoefficients<T, Basis> > save;

        if (residual<=tol && k>1) {
            #ifdef VERBOSE
                std::cout << "pcg: Tolerance reached r = "
                          << residual << std::endl;
            #endif
            return k;
        }

        /* alpha_k */
        tgal0 = std::chrono::system_clock::now();
        if (k==1) ak = residual*residual;
        else      ak = dot(r, p);
        tgal1 = std::chrono::system_clock::now();
        dgal  = std::chrono::duration_cast<std::chrono::seconds>(tgal1-tgal0);
        std::cout << "pcg: (r,p) time: " << dgal.count() << std::endl;

        tgal0 = std::chrono::system_clock::now();
        pAp = energy(stiff, S, p, p, Lambda, save);
        tgal1 = std::chrono::system_clock::now();
        dgal  = std::chrono::duration_cast<std::chrono::seconds>(tgal1-tgal0);
        std::cout << "pcg: (p, Ap) time: " << dgal.count() << std::endl;
        ak /= pAp;

        /* Update x_k and r_k */
        trunc_acc = nrmp*delta*std::fabs(ak);
        tgal0 = std::chrono::system_clock::now();
        if (uzero && k==1) {
            x.tree() = ak*p.tree();
        } else {
            tmp = p;
            scal(ak, tmp);
            x.tree() = add_truncate(tmp.tree(), x.tree(), trunc_acc);
        }
        tgal1 = std::chrono::system_clock::now();
        dgal  = std::chrono::duration_cast<std::chrono::seconds>(tgal1-tgal0);
        std::cout << "pcg: update x time: " << dgal.count() << std::endl;


        /* Compute residual */
        T restol   = residual/2.;
        tgal0 = std::chrono::system_clock::now();
        r          = evallaplace(stiff, S, x, Lambda, Lambda, restol/2.);
        tgal1 = std::chrono::system_clock::now();
        dgal  = std::chrono::duration_cast<std::chrono::seconds>(tgal1-tgal0);
        std::cout << "pcg: Ax time: " << dgal.count() << std::endl;
        scal(-1., r);
        tgal0 = std::chrono::system_clock::now();
        r.tree()   = add_truncate(b.tree(), r.tree(), restol/2.);
        tgal1 = std::chrono::system_clock::now();
        dgal  = std::chrono::duration_cast<std::chrono::seconds>(tgal1-tgal0);
        std::cout << "pcg: b-Ax time: " << dgal.count() << std::endl;
        residual   = nrm2(r);

        /* Update p_k */
        tgal0 = std::chrono::system_clock::now();
        rAp = energy(r, save);
        tgal1 = std::chrono::system_clock::now();
        dgal  = std::chrono::duration_cast<std::chrono::seconds>(tgal1-tgal0);
        std::cout << "pcg: (r, Ap) time: " << dgal.count() << std::endl;
        bk  = -rAp/pAp;
        tgal0 = std::chrono::system_clock::now();
        if (bk<=0) {
            p = r;
        } else {
            p.tree()  = r.tree()+bk*p.tree();
            nrmp      = nrm2(p);
            p.truncate(delta*nrmp);
        }
        tgal1 = std::chrono::system_clock::now();
        dgal  = std::chrono::duration_cast<std::chrono::seconds>(tgal1-tgal0);
        std::cout << "pcg: Update p time: " << dgal.count() << std::endl;

        #ifdef VERBOSE
            std::cout << "pcg: Iteration " << k
                      << " residual " << residual << std::endl;
            std::cout << "pcg: max rank r "
                      << r.tree().max_rank() << std::endl;
            std::cout << "pcg: max rank solution "
                      << x.tree().max_rank() << std::endl;
            std::cout << "pcg: max rank p "
                      << p.tree().max_rank() << std::endl;
            std::cout << "pcg: max rank b "
                      << b.tree().max_rank() << std::endl;
            std::cout << "pcg: trunc_acc    = " << trunc_acc
                      << std::endl;
            std::cout << "pcg: nrmp         = " << nrmp
                      << std::endl;
            std::cout << "pcg: ak           = " << ak
                      << std::endl;
            std::cout << "pcg: bk           = " << bk
                      << std::endl;
        #endif
    }

    std::cerr << "pcg: Max iterations reached: maxit " << maxit
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
                SeparableRHSD<T, Basis>& fint,
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
    std::vector<IndexSet<Index1D> > eval_diff(sweep.size());

    /* Determine extended index set for evaluation */
    for (size_type j=0; j<current.size(); ++j) {
        IndexSet<Index1D> currenteval = total[j];
        extendMultiTree(u.basis(), sweep[j], total[j],
                        "standard", false);
        diff[j]      = total[j];
        eval_diff[j] = total[j];

        for (auto& lambda : current[j]) {
            diff[j].erase(lambda);
        }

        for (auto& lambda : currenteval) {
            eval_diff[j].erase(lambda);
        }
    }

    /* Evaluate */
    genAddCoefficients(fcp, fint, eval_diff);
    set(f, fcp, total);

    /* Compute residual */
    r = eval(A, u, total, current);
    scal(-1., r);
    r.tree() = add_truncate(f.tree(), r.tree(), trunc);
    r        = applyScale(S, r, total, trunc);

    return diff;
}


template <typename Optype, typename Basis, typename T>
std::vector<IndexSet<Index1D> >
presidual2(      Sepop<Optype>& A,
                 Sepdiagscal<Basis>& S,
                 HTCoefficients<T, Basis>& u,
                 HTCoefficients<T, Basis>& f,
                 SepCoefficients<Lexicographical, T, Index1D>& fcp,
                 HTCoefficients<T, Basis>& r,
                 SeparableRHSD<T, Basis>& fint,
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
    std::vector<IndexSet<Index1D> > eval_diff(sweep.size());

    /* Determine extended index set for evaluation */
    auto start = std::chrono::system_clock::now();
    for (size_type j=0; j<current.size(); ++j) {
        IndexSet<Index1D> currenteval = total[j];
        extendMultiTree(u.basis(), sweep[j], total[j],
                        "standard", false);
        diff[j]      = total[j];
        eval_diff[j] = total[j];

        for (auto& lambda : current[j]) {
            diff[j].erase(lambda);
        }

        for (auto& lambda : currenteval) {
            eval_diff[j].erase(lambda);
        }
    }
    auto end     = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start);
    std::cout << "presidual: extend took " << elapsed.count() << " secs\n";

    /* Evaluate */
    start = std::chrono::system_clock::now();
    genAddCoefficients(fcp, fint, eval_diff);
    set(f, fcp, total);
    end     = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start);
    std::cout << "presidual: rhs took " << elapsed.count() << " secs\n";

    /* Compute residual */
    start = std::chrono::system_clock::now();
    r = evallaplace(A, S, u, total, current, trunc/3.);
    end     = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start);
    std::cout << "presidual: A*x took " << elapsed.count() << " secs\n";

    scal(-1., r);
    start = std::chrono::system_clock::now();
    f = applyScale(S, f, total, trunc/3.);
    end     = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start);
    std::cout << "presidual: S*f took " << elapsed.count() << " secs\n";

    start = std::chrono::system_clock::now();
    r.tree() = add_truncate(f.tree(), r.tree(), trunc/3.);
    end     = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start);
    std::cout << "presidual: Sf-Ax took " << elapsed.count() << " secs\n";

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
    size_type num = 0;
    T P_Lambda    = buckets.bucketsort2(contvec, thresh);

    thresh *= thresh;

    for (int i=buckets.bucket_ell2norms.size()-1; i>=0;
         --i) {
        P_Lambda += std::pow(buckets.bucket_ell2norms[i], 2.0L);
        if (P_Lambda > thresh) {
            num = i;
            P_Lambda -= std::pow(buckets.bucket_ell2norms[i], 2.0L);
            break;
        }
    }


    /* Add buckets */
    for (size_type i=0; i<=num; ++i) {
        buckets.addBucketToCoefficients(newind, i);
    }

    /* Remove indices from last bucket */
    Coefficients<Lexicographical, T, Index1DC>  remove;
    buckets.addBucketToCoefficients(remove, num);
    int count = 1;
    for (auto& lambda : remove) {
        P_Lambda += lambda.second*lambda.second;
        if (P_Lambda>thresh) {
            break;
        }
        if (count>1) newind.erase(lambda.first);
        ++count;
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
        auto resnew = nrm2(res);
        std::cout << "bulk: true alpha = " << resnew/resex << std::endl;
    #endif

    return sweep;
}


template <typename T, typename Basis>
std::vector<IndexSet<Index1D> >
bulknew(const T alpha, const T resex,
              HTCoefficients<T, Basis>& res,
              std::vector<IndexSet<Index1D> >& Lambda,
        const std::vector<IndexSet<Index1D> >& diff,
        const unsigned                         init_size = 1)
{
    assert(Lambda.size()==(unsigned) res.dim());

    typedef typename std::vector<IndexSet<Index1D> >::size_type size_type;

    std::vector<flens::DenseVector<flens::Array<T> > >
                                    sigmas(Lambda.size());
    std::vector<Coefficients<Lexicographical, T, Index1D> >
                                    cont(Lambda.size());
    std::vector<IndexSet<Index1D> > sweep(Lambda.size());

    T thresh    = alpha*resex;

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
    Coefficients<Lexicographical, T, Index1DC>  add;
    size_type scale     = 2;
    size_type size      = init_size;
    size_type count     = 0;
    size_type count_old = 0;
    T P_Lambda          = 0;
    buckets.bucketsort(contvec, std::sqrt(1.-alpha*alpha)*resex);

    size_type it = 1;
    for (unsigned i=0; i<buckets.bucket_ell2norms.size(); ++i) {
        /* Count number of indices */
        count_old  = count;
        count     += buckets.buckets[i].size();
        bool end = false;
        if (count<=size) {
            buckets.addBucketToCoefficients(add, i);
            if (i<buckets.buckets.size()-1) {
                continue;
            } else {
                end = true;
            }
        }

        count = count_old;
        for (const auto& lambda : buckets.buckets[i]) {
            if (count==size || end) {
                break;
            }

            ++count;
            add[(*lambda).first] = (*lambda).second;
        }

        /* Add indices */
        for (const auto& it : add) {
            completeMultiTree(res.basis(), it.first,
                              Lambda[it.first.d-1],
                              sweep[it.first.d-1],
                              true);
        }

        /* Compute restricted residual */
        auto copy = res;
        restrict(copy, Lambda);
        P_Lambda  = nrm2(copy);

        #ifdef VERBOSE
            std::cout << "bulk: Iteration " << it << ", true alpha = "
                      << P_Lambda/resex << std::endl;
        #endif

        if (P_Lambda >= thresh) {
            break;
        }

        /* Increase size */
        add.clear();
        ++it;
        size *= scale;
        count = count_old;
        if (!end) --i;
    }

    return sweep;
}


template <typename T, typename Basis>
std::vector<IndexSet<Index1D> >
bulk2(const T alpha, const T resex,
            HTCoefficients<T, Basis>& res,
            std::vector<IndexSet<Index1D> >& Lambda,
      const std::vector<IndexSet<Index1D> >& diff,
      const std::vector<IndexSet<Index1D> >& total)
{
    assert(Lambda.size()==(unsigned) res.dim());

    typedef typename std::vector<IndexSet<Index1D> >::size_type size_type;
    typedef typename flens::DenseVector<flens::Array<T> >       Vector;
    typedef Coefficients<Lexicographical, T, Index1D>           Coeff1D;

    std::vector<Vector>             sigmas(Lambda.size());
    std::vector<Coeff1D>            cont(Lambda.size());
    std::vector<IndexSet<Index1D> > sweep(Lambda.size());
    Vector                          base(Lambda.size());

    T thresh = std::sqrt((T) 1-alpha*alpha)*resex;

    // Compute contractions
    res.orthogonalize_svd(sigmas);
    for (int j=1; j<=base.length(); ++j) {
        auto copy   = res;
        auto active = total;
        active[j-1] = Lambda[j-1];
        restrict(copy, active);
        T resbase   = nrm2(copy);
        base(j)     = resbase;
        std::cout << "Minimum = " << resbase/resex << std::endl;
    }
    contraction(res, diff, sigmas, cont);


    // Compute index set
    T stop = alpha*alpha*resex*resex;
    for (size_type j=0; j<Lambda.size(); ++j) {
        // Bucket sort
        Coefficients<Bucket, T, Index1D>          buckets;
        Coefficients<Lexicographical, T, Index1D> newind;
        buckets.bucketsort(cont[j], thresh);

        // Determine bulk
        T P_Lambda = base(j+1)*base(j+1);
        size_type i = 0;
        std::cout << "res/res = "
                  << std::sqrt(base(j+1)*base(j+1)
                  +cont[j].norm(2.)*cont[j].norm(2.))/resex << std::endl;
        for (; i<buckets.bucket_ell2norms.size()-1; ++i) {
            P_Lambda += std::pow(buckets.bucket_ell2norms[i], 2.0L);
            if (P_Lambda >= stop) {
                P_Lambda -= std::pow(buckets.bucket_ell2norms[i], 2.0L);
                break;
            }
            (void) buckets.addBucketToCoefficients(newind, i);
        }

        for (const auto& it : buckets.buckets[i]) {
            if (P_Lambda >= stop) break;
            newind[it->first] = it->second;
            P_Lambda         += it->second*it->second;
        }

        std::cout << "newind/res = "
                  << std::sqrt(base(j+1)*base(j+1)
                  +newind.norm(2.)*newind.norm(2.))/resex << std::endl;

        // Complete tree
        for (const auto& it : newind) {
            completeMultiTree(res.basis(), it.first,
                              Lambda[j],
                              sweep[j], true);
        }
    }

    #ifdef VERBOSE
        restrict(res, Lambda);
        auto resnew = nrm2(res);
        std::cout << "bulk: true alpha = " << resnew/resex << std::endl;
    #endif

    return sweep;
}


template <typename T, typename Basis>
std::vector<IndexSet<Index1D> >
bulkBestN(const T alpha, const T resex,
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

    /* Sort */
    Coefficients<AbsoluteValue, T, Index1DC>    sorted;
    Coefficients<Lexicographical, T, Index1DC>  newind;
    sorted = contvec;
    T P_Lambda    = 0;

    thresh *= thresh;

    auto it = sorted.end();
    --it;
    auto save = it;
    for (;it!=sorted.begin(); --it) {
        P_Lambda += std::pow((*it).first, 2.0L);
        if (P_Lambda > thresh) {
            break;
        }
        --save;
    }

    /* Add buckets */
    for (it=sorted.begin();; ++it) {
        newind[(*it).second] = (*it).first;
        if (it==save) break;
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
        auto resnew = nrm2(res);
        std::cout << "bulk: true alpha = " << resnew/resex << std::endl;
    #endif

    return sweep;
}


template <typename Optype, typename Basis, typename T>
unsigned
htawgm(Sepop<Optype>&                   A,
       Sepdiagscal<Basis>&              S,
       HTCoefficients<T, Basis>&        u,
       SeparableRHSD<T, Basis>&   f,
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

    /* Initial residual */
    auto start = std::chrono::system_clock::now();
    if (!params.uzero) {
        Au = eval(A, u, Lambda, Lambda);
        Au.truncate(nrmf*1e-01);
        scal(-1., Au);
        r.tree() = add_truncate(F.tree(), Au.tree(), nrmf*1e-01);
    }

    r         = applyScale(S, r, Lambda, nrmf*1e-01);
    residual  = nrm2(r);

    if (residual<=params.tol_awgm) {
        #ifdef VERBOSE
            std::cout << "htawgm: Tolerance reached r = " << residual
                      << std::endl;
        #endif

        return 0;
    }

    T tol;
    T gamma = params.gamma0;
    for (unsigned k=1; k<=params.maxit_awgm; ++k) {
        unsigned pcg_it, size;
        T        res_pcg;

        auto end     = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start);
        std::cout << "htawgm: Iteration " << k
                  << ", current elapsed time "
                  << elapsed.count() << " seconds\n";

        #ifdef VERBOSE
            std::cout << "htawgm: Iteration " << k
                      << " r = " << residual << std::endl;
        #endif

        /* Galerkin solve */
        tol    = gamma*residual;
        pcg_it = galerkin_pcg(A, S, u, F, Lambda, res_pcg,
                               params.uzero,
                               tol,
                               params.maxit_pcg,
                               params.delta1_pcg,
                               params.delta2_pcg,
                               params.delta3_pcg,
                               1e-02);
        #ifdef VERBOSE
            std::cout << "htawgm: galerkin_pcg required " << pcg_it
                      << " iterations to reach tolerance "
                      << tol
                      << std::endl;
        #else
            (void) pcg_it;
        #endif

        /* Approximate residual */
        T save       = S.nu();
        T cv         = (1.-S.eps())/(4.*params.nrmA);
        T eta        = 0.5*cv*params.omega*residual/nrm2(F);
        S.set_nu(eta);
        sweep        = presidual(A, S, u, F, Fcp, r, f,
                                 Lambda, sweep, total,
                                 params.tol_awgm);
        S.set_nu(save);
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
            std::cout << "htawgm: max rank solution " << u.tree().max_rank()
                      << std::endl;
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


template <typename Optype, typename Basis, typename T>
unsigned
richardson(Sepop<Optype>&                   A,
           Sepdiagscal<Basis>&              S,
           HTCoefficients<T, Basis>&        x,
           HTCoefficients<T, Basis>&        b,
           std::vector<IndexSet<Index1D> >& Lambda,
           T                                omega,
           T                                tol,
           unsigned                         maxit)
{
    auto r     = b;
    T residual = nrm2(r);
    T delta    = 1e-01;

    std::cout << "richardson: Iteration 0, residual "
              << residual << std::endl;
    if (residual<=tol) return 0;
    auto stiff = assemble_laplace<T, Optype>(A, Lambda, Lambda);

    for (unsigned k=1; k<=maxit; ++k) {
        scal(omega, r);
        T trunc  = delta*omega*residual;
        if (k==1) x   = r;
        else x.tree() = add_truncate(x.tree(), r.tree(), trunc);

        r        = evallaplace(stiff, S, x, Lambda, Lambda, trunc/2.);
        scal(-1., r);
        r.tree() = add_truncate(b.tree(), r.tree(), trunc/2.);
        residual = nrm2(r);

        std::cout << "richardson: Iteration " << k
                  << ", residual " << residual << std::endl;
        std::cout << "richardson: rank residual " << r.tree().max_rank()
                  << std::endl;
        std::cout << "richardson: rank solution " << x.tree().max_rank()
                  << std::endl;
        if (residual<=tol) return k;
    }

    return maxit;
}


template <typename Optype, typename Basis, typename T>
unsigned
htawgm2(      Sepop<Optype>&                   A,
              Sepdiagscal<Basis>&              S,
              HTCoefficients<T, Basis>&        u,
              SeparableRHSD<T, Basis>&   f,
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
    HTCoefficients<T, Basis>                     SF(u.dim(),    u.basis(),
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

    /* Initial residual */
    r  = applyScale(S, F, Lambda, params.tol_awgm);
    SF = r;
    auto start = std::chrono::system_clock::now();
    if (!params.uzero) {
        Au = evaleff(A, S, u, Lambda, Lambda, params.tol_awgm);
        scal(-1., Au);
        r.tree() = add_truncate(r.tree(), Au.tree(), params.tol_awgm);
    }

    residual  = nrm2(r);
    #ifdef VERBOSE
        std::cout << "htawgm: Iteration " << 0
                  << " r = " << residual << std::endl;
    #endif


    if (residual<=params.tol_awgm) {
        #ifdef VERBOSE
            std::cout << "htawgm: Tolerance reached r = " << residual
                      << std::endl;
        #endif

        return 0;
    }

    T tol;
    T gamma   = params.gamma0;
    T ksi     = residual;
    T omega3  = 1e-02;
    T omega4  = 1.;
    for (unsigned k=1; k<=params.maxit_awgm; ++k) {
        unsigned pcg_it, size;
        T        res_pcg;

        auto end     = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start);

        /* Galerkin solve */
        tol    = gamma*residual;
        auto tgal0 = std::chrono::system_clock::now();
        pcg_it = galerkin_pcg2(A, S, u, SF, Lambda, res_pcg,
                               params.uzero,
                               tol,
                               params.maxit_pcg,
                               params.delta1_pcg,
                               residual/2.);
        auto tgal1 = std::chrono::system_clock::now();
        auto dgal  = std::chrono::duration_cast<std::chrono::seconds>(tgal1-tgal0);
        std::cout << "galerkin time: " << dgal.count() << std::endl;
        #ifdef VERBOSE
            std::cout << "htawgm: galerkin_pcg required " << pcg_it
                      << " iterations to reach tolerance "
                      << tol
                      << std::endl;
        #else
            (void) pcg_it;
        #endif

        /* Approximate residual */
        tgal0 = std::chrono::system_clock::now();
        sweep      = presidual2(A, S, u, F, Fcp, r, f,
                                Lambda, sweep, total,
                                residual*gamma);
        tgal1 = std::chrono::system_clock::now();
        dgal  = std::chrono::duration_cast<std::chrono::seconds>(tgal1-tgal0);
        std::cout << "presidual time: " << dgal.count() << std::endl;

        tgal0 = std::chrono::system_clock::now();
        residual = nrm2(r);
        tgal1 = std::chrono::system_clock::now();
        dgal  = std::chrono::duration_cast<std::chrono::seconds>(tgal1-tgal0);
        std::cout << "residual norm time: " << dgal.count() << std::endl;
        std::cout << "htawgm: Iteration " << k
                  << ", current elapsed time "
                  << elapsed.count() << " seconds\n";

        #ifdef VERBOSE
            std::cout << "htawgm: Iteration " << k
                      << " r = " << residual << std::endl;
        #endif

        if (residual<=params.tol_awgm) {
            #ifdef VERBOSE
                std::cout << "htawgm: Tolerance reached r = " << residual
                          << std::endl;
            #endif

            return k;
        }

        /* Recompress */
        if (residual<=omega3*ksi) {
            u.truncate(omega4*residual);
            T trunc = residual*gamma;
            r       = evallaplace(A, S, u, total, Lambda, trunc);
            scal(-1., r);
            auto tmp = applyScale(S, F, total, trunc);
            r.tree() = add_truncate(tmp.tree(), r.tree(), trunc);
            residual = nrm2(r);
            ksi      = residual;
            #ifdef VERBOSE
                std::cout << "htawgm: Iteration " << k
                          << ", recompress, r = " << residual << std::endl;
            #endif
        }

        /* Bulk chasing */
        tgal0 = std::chrono::system_clock::now();
        sweep = bulknew(params.alpha, residual,
                        r, Lambda, sweep);
        tgal1 = std::chrono::system_clock::now();
        dgal  = std::chrono::duration_cast<std::chrono::seconds>(tgal1-tgal0);
        std::cout << "bulk time: " << dgal.count() << std::endl;

        /* New RHS */
        restrict(F, Lambda);
        SF = applyScale(S, F, Lambda, tol);

        /* Extend u to new Lambda */
        extend(u, Lambda);
        params.uzero = false;

        #ifdef VERBOSE
            std::cout << "htawgm: max rank solution " << u.tree().max_rank()
                      << std::endl;
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
