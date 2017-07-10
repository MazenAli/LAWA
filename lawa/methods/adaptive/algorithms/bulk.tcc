#ifndef LAWA_METHODS_ADAPTIVE_ALGORITHMS_BULK_TCC
#define LAWA_METHODS_ADAPTIVE_ALGORITHMS_BULK_TCC 1

#include <cmath>
#include <cassert>
#include <iostream>

#include <lawa/methods/adaptive/algorithms/multitreeoperations.h>

namespace lawa
{

template <typename T, typename Basis>
IndexSet<Index1D>
bulk(      IndexSet<Index1D>&                         active,
           Coefficients<Lexicographical, T, Index1D>& r,
     const T                                          alpha,
     const T                                          res,
     const T                                          res_base,
     const Basis&                                     basis,
     const bool                                       verbose)
{
    assert(alpha>0. && alpha<=1.);

    // Set up difference vector
    Coefficients<Lexicographical, T, Index1D> rdiff;
    for (const auto& it : r)
        if (active.find(it.first)==active.end()) rdiff[it.first] = it.second;

    // Parameters
    T thresh = std::sqrt(1.-alpha*alpha)*res;
    T stop   = alpha*alpha*res*res;

    // Bucket sort
    Coefficients<Bucket, T, Index1D>          buckets;
    Coefficients<Lexicographical, T, Index1D> newind;
    buckets.bucketsort(rdiff, thresh);
    T P_Lambda = res_base*res_base;

    // Add coefficients
    unsigned i=0;
    for (; i<buckets.bucket_ell2norms.size()-1; ++i) {
        P_Lambda += std::pow(buckets.bucket_ell2norms[i], 2.);
        if (P_Lambda >= stop) {
            P_Lambda -= std::pow(buckets.bucket_ell2norms[i], 2.);
            break;
        }
        (void) buckets.addBucketToCoefficients(newind, i);
    }

    for (const auto& it : buckets.buckets[i]) {
        if (P_Lambda >= stop) break;
        newind[it->first] = it->second;
        P_Lambda         += it->second*it->second;
    }

    // Complete tree
    IndexSet<Index1D> sweep;
    for (const auto& it : newind) {
        completeMultiTree(basis, it.first, active, sweep, true);
    }

    if (verbose) {
        Coefficients<Lexicographical, T, Index1D> rres;
        for (const auto& it : active) {
            rres[it] = r[it];
        }

        T a = rres.norm(2.)/res;
        std::cout << "bulk: alpha = " << a << std::endl;
    }

    return sweep;
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_BULK_TCC
