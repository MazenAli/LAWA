#ifndef LAWA_METHODS_ADAPTIVE_PRECONDITIONERS_DIAGONALPRECONDITIONER_TCC
#define LAWA_METHODS_ADAPTIVE_PRECONDITIONERS_DIAGONALPRECONDITIONER_TCC 1

namespace lawa
{

template <typename P, typename T>
DiagonalPreconditioner<P, T>::
DiagonalPreconditioner(P& _p):
    p_(_p)
{}


template <typename P, typename T>
void
DiagonalPreconditioner<P, T>::
eval(Coefficients<Lexicographical, T, Index1D>& v)
{
    for (auto& it : v) {
        it.second *= p_(it.first);
    }
}


template <typename P, typename T>
void
DiagonalPreconditioner<P, T>::
operator()(Coefficients<Lexicographical, T, Index1D>& v)
{
    return eval(v);
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_PRECONDITIONERS_DIAGONALPRECONDITIONER_TCC
