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
         const std::vector<IndexSet<Index1D> >& Lambda,
               std::vector<IndexSet<Index1D> >& sweep,
               std::vector<IndexSet<Index1D> >& total)
{
    assert(A.dim()==(unsigned) u.dim());
    assert(A.dim()==(unsigned) f.dim());
    assert(A.dim()==Lambda.size());
    assert(A.dim()==sweep.size());
    assert(A.dim()==total.size());

    typedef typename std::vector<IndexSet<Index1D> >::size_type size_type;
    std::vector<IndexSet<Index1D> > eval_diff(sweep.size());

    /* Determine extended index set for evaluation */
    for (size_type j=0; j<Lambda.size(); ++j) {
        IndexSet<Index1D> currenteval = total[j];
        extendMultiTree(u.basis(), sweep[j], total[j],
                        "standard", false);
        sweep[j]     = total[j];
        eval_diff[j] = total[j];

        for (auto& lambda : Lambda[j]) {
            sweep[j].erase(lambda);
        }

        for (auto& lambda : currenteval) {
            eval_diff[j].erase(lambda);
        }
    }

    /* Evaluate */
    genAddCoefficients(fcp, fint, eval_diff);
    set(f, fcp);

    /* Compute residual */
    HTCoefficients<T, Basis> r = eval(A, u, total, Lambda);
    scal(-1., r);
    r.tree() = f.tree()+r.tree();

    return r;
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_SZONERES_TCC
