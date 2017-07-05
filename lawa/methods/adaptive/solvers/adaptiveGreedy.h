#ifndef LAWA_METHODS_ADAPTIVE_SOLVERS_ADAPTIVEGREEDY_H
#define LAWA_METHODS_ADAPTIVE_SOLVERS_ADAPTIVEGREEDY_H 1

namespace lawa
{

unsigned
adaptiveGreedy(Engine* ep,
                Sepop<Optype>&  A,
                Sepdiagscal<Basis>& S,
                Prec&               P,
                HTCoefficients<T, Basis> u,


} // namespace lawa

#include <lawa/methods/adaptive/solvers/adaptiveGreedy.tcc>

#endif // LAWA_METHODS_ADAPTIVE_SOLVERS_ADAPTIVEGREEDY_H
