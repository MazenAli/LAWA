namespace lawa {

template <typename T>
MRA<T,Orthogonal,R,MultiRefinement>::MRA(FLENS_DEFAULT_INDEXTYPE _d, FLENS_DEFAULT_INDEXTYPE j)
    : d(_d), j0(j), _j(j0), phi(d)
{
    if (d > 4) {
        std::cerr << "MRA<T,Orthogonal,R,MultiRefinement> not yet implemented for d = " << d << std::endl;
        exit(1);
    }
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Orthogonal,R,MultiRefinement>::level() const
{
    return _j;
}

template <typename T>
void
MRA<T,Orthogonal,R,MultiRefinement>::setLevel(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(d==1 || j>=j0);
    _j = j;
}

} // namespace lawa

