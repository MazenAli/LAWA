#ifndef LAWA_METHODS_ADAPTIVE_ALGORITHMS_SAMPLE_H
#define LAWA_METHODS_ADAPTIVE_ALGORITHMS_SAMPLE_H


#include <cstddef>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/settings/enum.h>

namespace lawa
{


template <typename T, typename _Basis, typename _Index, typename _Rhs>
void
sample_f(const _Basis& basis, IndexSet<_Index> Lambda, _Rhs& f,
         Coefficients<Lexicographical, T, _Index>& ret,
         T tol,
         bool IsMW = false, T alpha = 0.7, std::size_t max_it = 1e+02);


template <typename T, typename _Basis, typename _Index,
         typename _Rhs, typename _Precon>
void
sample_f(const _Basis& basis, IndexSet<_Index> Lambda, _Rhs& f,
         const _Precon& P,
         Coefficients<Lexicographical, T, _Index>& ret,
         T tol,
         bool IsMW = false, T alpha = 0.7, std::size_t max_it = 1e+02);


// u assumed to be properly scaled
template <typename T, typename _Basis, typename _Index, typename _Op,
          typename _Precon>
void
sample_Au(const _Basis& basis, IndexSet<_Index> Lambda, _Op& A,
          const _Precon& P,
          const Coefficients<Lexicographical, T, _Index>& u,
          Coefficients<Lexicographical, T, _Index>& ret,
          T tol,
          bool IsMW = false, T alpha = 0.7, std::size_t max_it = 1e+02);


} // namespace lawa


#include <lawa/methods/adaptive/algorithms/sample.tcc>

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_SAMPLE_H
