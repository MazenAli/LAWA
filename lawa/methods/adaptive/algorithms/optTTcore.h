#ifndef LAWA_METHODS_ADAPTIVE_ALGORITHMS_OPTTTCORE_H
#define LAWA_METHODS_ADAPTIVE_ALGORITHMS_OPTTTCORE_H 1

#include <vector>
#include <engine.h>
#include <flens/flens.cxx>

#include <lawa/methods/adaptive/solvers/solver_parameters.h>

#ifdef VERBOSE
    #ifndef BUFSIZE
        #define BUFSIZE 524288
    #endif
#endif

namespace lawa
{

/* Optimize HT-core of a TT tree using Anthony's Matlab toolbox.
 * Laplace-like operator.
 */
template <typename T, typename I>
std::vector<flens::GeMatrix<
flens::FullStorage<T, cxxblas::ColMajor> > >
optTTcoreLaplace(      Engine                                       *ep,
                 const std::vector<flens::GeMatrix<
                       flens::FullStorage<T, cxxblas::ColMajor> > >& A,
                 const std::vector<flens::GeMatrix<
                       flens::FullStorage<T, cxxblas::ColMajor> > >& rhs,
                 const std::vector<flens::GeMatrix<
                       flens::FullStorage<T, cxxblas::ColMajor> > >& x0,
                       flens::DenseVector<
                       flens::Array<I> >&                            ranks,
                 const OptTTCoreParams&                              params);

} // namespace lawa

#include <lawa/methods/adaptive/algorithms/optTTcore.tcc>

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_OPTTTCORE_H
