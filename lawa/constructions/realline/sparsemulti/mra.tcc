namespace lawa {

template <typename T>
MRA<T,Primal,R,SparseMulti>::MRA(FLENS_DEFAULT_INDEXTYPE _d, FLENS_DEFAULT_INDEXTYPE j)
    : d(_d), j0(j), _j(j0), phi(d)
{

}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Primal,R,SparseMulti>::level() const
{
    return _j;
}

template <typename T>
void
MRA<T,Primal,R,SparseMulti>::setLevel(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=j0);
    _j = j;
}

}   // namespace lawa
