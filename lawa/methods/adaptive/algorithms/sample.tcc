#ifndef LAWA_METHODS_ADAPTIVE_ALGORITHMS_SAMPLE_TCC
#define LAWA_METHODS_ADAPTIVE_ALGORITHMS_SAMPLE_TCC


#include <cstddef>
#include <cassert>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/algorithms/multitreeoperations.h>
#include <cmath>


namespace lawa
{


template <typename T, typename _Basis, typename _Index, typename _Rhs>
void
sample_f(const _Basis& basis, IndexSet<_Index> Lambda, _Rhs& f,
         Coefficients<Lexicographical, T, _Index>& ret,
         T tol,
         bool IsMW, T alpha, std::size_t max_it)
{

    assert(Lambda.size() > 0);

    typedef     IndexSet<_Index>                            IndexSet;
    typedef     Coefficients<Lexicographical, T, _Index>    _Coefficients;

    _Coefficients    sweep;  // current extension sweep direction
    _Coefficients    total;  // total extention residual

    // Initial step
    FillWithZeros(Lambda, sweep);
    FillWithZeros(Lambda, total);
    ret.clear();
    ret = f(Lambda);

    // Sampling loop
    for (std::size_t k = 0; k < max_it; ++k) {
        extendMultiTree(basis, sweep, total, "standard", IsMW, true);
        IndexSet diff = supp(total);

        for (auto& lambda : Lambda) {
            diff.erase(lambda);
        }

        _Coefficients res = f(diff);
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
        Coefficients<Bucket, T, _Index> buckets;
        buckets.bucketsort(res, thresh);

        _Coefficients new_ind;
        for (std::size_t i = 0; i < buckets.bucket_ell2norms.size(); ++i) {
            P_Lambda += std::pow(buckets.bucket_ell2norms[i], 2.0L);
            buckets.addBucketToCoefficients(new_ind, i);

            if (P_Lambda >= std::pow(alpha*res_norm, (T)2.)) break;
        }

        // New index set
        IndexSet new_sweep;

        for (auto& coeff : new_ind) {
            if (ret.find(coeff.first) == ret.end()) {
                completeMultiTree(basis, coeff.first, ret,
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
        FillWithZeros(new_sweep, sweep);
        ret += f(new_sweep);
        Lambda = supp(ret);
    }

    std::cerr << "Warning! Max iterations " << max_it << " reached!\n";
}


template <typename T, typename _Basis, typename _Index,
         typename _Rhs, typename _Precon>
void
sample_f(const _Basis& basis, IndexSet<_Index> Lambda, _Rhs& f,
         _Precon& P,
         Coefficients<Lexicographical, T, _Index>& ret,
         T tol,
         bool IsMW, T alpha, std::size_t max_it)
{

    assert(Lambda.size() > 0);

    typedef     IndexSet<_Index>                            IndexSet;
    typedef     Coefficients<Lexicographical, T, _Index>    _Coefficients;

    _Coefficients    sweep;  // current extension sweep direction
    _Coefficients    total;  // total extention residual

    // Initial step
    FillWithZeros(Lambda, sweep);
    FillWithZeros(Lambda, total);
    ret.clear();
    ret = f(Lambda);

    // Scale
    for (auto& lambda : ret) {
        lambda.second *= P(lambda.first);
    }

    // Sampling loop
    for (std::size_t k = 0; k < max_it; ++k) {
        extendMultiTree(basis, sweep, total, "standard", IsMW, true);
        IndexSet diff = supp(total);

        for (auto& lambda : Lambda) {
            diff.erase(lambda);
        }

        _Coefficients res = f(diff);

        for (auto& lambda : res) {
            lambda.second *= P(lambda.first);
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
        Coefficients<Bucket, T, _Index> buckets;
        buckets.bucketsort(res, thresh);

        _Coefficients new_ind;
        for (std::size_t i = 0; i < buckets.bucket_ell2norms.size(); ++i) {
            P_Lambda += std::pow(buckets.bucket_ell2norms[i], 2.0L);
            buckets.addBucketToCoefficients(new_ind, i);

            if (P_Lambda >= std::pow(alpha*res_norm, (T)2.)) break;
        }

        // New index set
        IndexSet new_sweep;

        for (auto& coeff : new_ind) {
            if (ret.find(coeff.first) == ret.end()) {
                completeMultiTree(basis, coeff.first, ret,
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
        FillWithZeros(new_sweep, sweep);
        ret += f(new_sweep);

        for (auto& lambda : new_sweep) {
            ret[lambda] *= P(lambda);
        }

        Lambda = supp(ret);
    }

    std::cerr << "Warning! Max iterations " << max_it << " reached!\n";
}


template <typename T, typename _Basis, typename _Index, typename _Op,
          typename _Precon>
void
sample_Au(const _Basis& basis, IndexSet<_Index> Lambda, _Op& A,
          _Precon& P,
          const Coefficients<Lexicographical, T, _Index>& u,
          Coefficients<Lexicographical, T, _Index>& ret, T tol,
          bool IsMW, T alpha, std::size_t max_it)
{

    assert(Lambda.size() > 0);

    typedef     IndexSet<_Index>                            IndexSet;
    typedef     Coefficients<Lexicographical, T, _Index>    _Coefficients;

    _Coefficients    sweep;  // current extension sweep direction
    _Coefficients    total;  // total extention residual

    // Initial step
    FillWithZeros(Lambda, sweep);
    FillWithZeros(Lambda, total);
    ret.clear();
    FillWithZeros(Lambda, ret);
    A.eval(u, ret);

    // Scale
    for (auto& lambda : ret) {
        lambda.second *= P(lambda.first);
    }

    // Sampling loop
    for (std::size_t k = 0; k < max_it; ++k) {
        extendMultiTree(basis, sweep, total, "standard", IsMW, true);
        _Coefficients res = total;
        A.eval(u, res);

        for (auto& lambda : res) {
            lambda.second *= P(lambda.first);
        }

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
        Coefficients<Bucket, T, _Index> buckets;
        buckets.bucketsort(res, thresh);

        _Coefficients new_ind;
        for (std::size_t i = 0; i < buckets.bucket_ell2norms.size(); ++i) {
            P_Lambda += std::pow(buckets.bucket_ell2norms[i], 2.0L);
            buckets.addBucketToCoefficients(new_ind, i);

            if (P_Lambda >= std::pow(alpha*res_norm, (T)2.)) break;
        }

        // New index set
        IndexSet new_sweep;

        for (auto& coeff : new_ind) {
            if (ret.find(coeff.first) == ret.end()) {
                completeMultiTree(basis, coeff.first, ret,
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
        FillWithZeros(new_sweep, sweep);
        ret.setToZero();
        A.eval(u, ret);

        for (auto& lambda : ret) {
            lambda.second *= P(lambda.first);
        }

        Lambda = supp(ret);
    } // sampling loop

    std::cerr << "Warning! Max iterations " << max_it << " reached!\n";
}


} // namespace lawa


#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_SAMPLE_TCC
