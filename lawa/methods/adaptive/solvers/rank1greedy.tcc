#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_RANK1GREEDY_TCC
#define LAWA_METHODS_ADAPTIVE_SOLVERS_RANK1GREEDY_TCC 1

#include <cassert>

namespace lawa
{

/* Rank 1 one greedy (A s.p.d.) over fixed wavelet index set
 * Input current rank
 */
template <typename Optype, typename Prec, typename T, typename Basis>
void
rank1greedy_sym(      Sepop<Optype>&                      A,
                      Prec&                               P,
                      HTCoefficients<T, Basis>&           x,
                const HTCoefficients<T, Basis>&           b,
                const std::vector<IndexSet<Index1D> >&    Lambda,
                      T&                                  residual,
                const unsigned                            numUps,
                const Rank1UP_Params&                     params)
{
    assert(A.dim()==(unsigned)x.dim());
    assert(A.dim()==(unsigned)b.dim());
    assert(A.dim()==Lambda.size());

    typedef HTCoefficients<T, Basis>            HTCoeff;

    Rank1UP_Params  p = params;
    unsigned nsweeps;

    auto nrmb = nrm2(const_cast<HTCoeff&>(b));
    for (unsigned k=1; k<=numUps; ++k) {
        nsweeps = rank1update_sym(A, P, x, b, Lambda, p);
        #ifdef VERBOSE
            HTCoeff Ax = eval(A, x, Lambda, Lambda);
            Ax.tree()  = b.tree()-Ax.tree();
            residual   = nrm2(Ax)/nrmb;
            std::cout << "rank1greedy_sym: On update " << k
                      << " ALS required " << nsweeps
                      << " sweeps to reach relative residual "
                      << residual << std::endl;
        #endif
        p.update = true;
    }
    #ifndef VERBOSE
        HTCoeff Ax = eval(A, x, Lambda, Lambda);
        Ax.tree()  = b.tree()-Ax.tree();
        residual   = nrm2(Ax)/nrmb;
    #endif
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_RANK1GREEDY_TCC
