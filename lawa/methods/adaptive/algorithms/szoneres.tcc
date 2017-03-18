#ifndef LAWA_METHODS_ADAPTIVE_ALGORITHMS_SZONERES_TCC
#define LAWA_METHODS_ADAPTIVE_ALGORITHMS_SZONERES_TCC 1

namespace lawa
{

template <typename Optype, typename Basis, typename T>
HTCoefficients<T, Basis>
szoneres(      Sepop<Optype>& A,
               HTCoefficients<T, Basis>& u,
               HTCoefficients<T, Basis>& f,
               SepCoefficients<Lexicographical, T, Index1D>& fcp,
         const SeparableRHSD<T, Basis>& fint,
         const std::vector<IndexSet<Index1D> >& current,
               std::vector<IndexSet<Index1D> >& sweep,
               std::vector<IndexSet<Index1D> >& total)
{
    assert(A.dim()==(unsigned) u.dim());
    assert(A.dim()==(unsigned) f.dim());
    assert(A.dim()==current.size());
    assert(A.dim()==sweep.size());
    assert(A.dim()==total.size());

    typedef typename std::vector<IndexSet<Index1D> >::size_type size_type;

    /* Determine extended index set for evaluation */
    for (size_type j=0; j<current.size(); ++j) {
        extendMultiTree(u.basis(), sweep[j], total[j],
                        "standard", false);
        sweep[j] = total[j];

        for (auto& lambda : current[j]) {
            sweep[j].erase(lambda);
        }
    }

    /* Evaluate */
    genAddCoefficients(fcp, fint, sweep);
    set(f, fcp);

    /* Compute residual */
    HTCoefficients<T, Basis> r = eval(A, u, total, current);
    scal(-1., r);
    r.tree() = f.tree()+r.tree();

    return r;
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_SZONERES_TCC
