#ifndef LAWA_METHODS_ADAPTIVE_PRECONDITIONERS_DIAGONALPRECONDITIONER_H
#define LAWA_METHODS_ADAPTIVE_PRECONDITIONERS_DIAGONALPRECONDITIONER_H 1

#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/datastructures/index.h>

namespace lawa
{

// A wrapper for index-wise scaling
template <typename P, typename T>
class DiagonalPreconditioner
{

public:
    DiagonalPreconditioner(P& _p);

    void
    eval(Coefficients<Lexicographical, T, Index1D>& v);

    void
    operator()(Coefficients<Lexicographical, T, Index1D>& v);

private:
    P& p_;

};

} // namespace lawa

#include <lawa/methods/adaptive/preconditioners/diagonalPreconditioner.tcc>

#endif // LAWA_METHODS_ADAPTIVE_PRECONDITIONERS_DIAGONALPRECONDITIONER_H
