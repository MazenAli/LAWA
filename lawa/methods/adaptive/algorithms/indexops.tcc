#ifndef LAWA_METHODS_ADAPTIVE_ALGORITHMS_INDEXOPS_TCC
#define LAWA_METHODS_ADAPTIVE_ALGORITHMS_INDEXOPS_TCC 1

#include <iostream>
#include <stdlib.h>
#include <cassert>
#include <lawa/settings/enum.h>
#include <limits>

namespace lawa
{

template <typename Index, typename Basis>
unsigned FLENS_DEFAULT_INDEXTYPE
maptoint(const Index&, const Basis&)
{
    std::cerr << "error: maptoint not implemented for given Index type\n";
    exit(EXIT_FAILURE);
}


template <typename Basis>
unsigned FLENS_DEFAULT_INDEXTYPE
maptoint(const Index1D& index, const Basis& basis)
{
    FLENS_DEFAULT_INDEXTYPE j       = index.j;
    FLENS_DEFAULT_INDEXTYPE k       = index.k;
    XType type   = index.xtype;
    FLENS_DEFAULT_INDEXTYPE j0      = basis.j0;
    FLENS_DEFAULT_INDEXTYPE offset  = 0;

    assert(type==XBSpline || type==XWavelet);
    if (type == XBSpline) {
        assert(j==j0);
        assert(k>=basis.mra.rangeI(j).firstIndex() &&
               k<=basis.mra.rangeI(j).lastIndex());
        offset -= basis.mra.rangeI(j).firstIndex();
        ++offset;
        return k + offset;
    } else {
        assert(j>=j0);
        assert(k>=basis.rangeJ(j).firstIndex() &&
               k<=basis.rangeJ(j).lastIndex());
        offset += basis.mra.cardI(j);
        offset -= basis.rangeJ(j).firstIndex();
        ++offset;
        return k + offset;
    }
}


template <typename Index, typename Basis>
unsigned FLENS_DEFAULT_INDEXTYPE
maxintind(const IndexSet<Index>& indexset, const Basis& basis)
{
    unsigned FLENS_DEFAULT_INDEXTYPE max = 0;
    for (const auto& it : indexset) {
        unsigned FLENS_DEFAULT_INDEXTYPE idx = maptoint(it, basis);
        max = (idx>max) ? idx : max;
    }

    return max;
}


template <typename Index, typename T, typename Basis>
unsigned FLENS_DEFAULT_INDEXTYPE
maxintindhash(const IndexSet<Index>& indexset,
              const int dim,
              HTCoefficients<T, Basis>& u)
{
    unsigned FLENS_DEFAULT_INDEXTYPE max = 0;
    for (const auto& it : indexset) {
        unsigned FLENS_DEFAULT_INDEXTYPE idx = u.map()(it, dim);
        max = (idx>max) ? idx : max;
    }

    return max;
}


template <typename Index, typename Basis>
Index
maptowav(const unsigned FLENS_DEFAULT_INDEXTYPE, const Basis&)
{
    std::cerr << "error: maptowav not implemented for given Index type\n";
    exit(EXIT_FAILURE);
}


template <typename Basis>
Index1D
maptowav(const unsigned FLENS_DEFAULT_INDEXTYPE i, const Basis& basis)
{
    assert(i>0);

    FLENS_DEFAULT_INDEXTYPE j0 = basis.j0;
    FLENS_DEFAULT_INDEXTYPE j, k;
    XType type;
    unsigned FLENS_DEFAULT_INDEXTYPE bound;

    bound = basis.mra.rangeI(j0).lastIndex()-
            basis.mra.rangeI(j0).firstIndex()+1;
    if (i<=bound) {
        type = XBSpline;
        j = j0;
        k = i+basis.mra.rangeI(j).firstIndex()-1;
        return Index1D(j, k, type);
    }

    type = XWavelet;
    for (FLENS_DEFAULT_INDEXTYPE _j=j0; _j<=std::numeric_limits<FLENS_DEFAULT_INDEXTYPE>::max(); ++_j) {
        assert(bound==(unsigned) basis.mra.cardI(_j));
        bound += basis.rangeJ(_j).lastIndex()-
                 basis.rangeJ(_j).firstIndex()+1;
        if (i<=bound) {
            j = _j;
            k = i+basis.rangeJ(j).firstIndex()-basis.mra.cardI(j)-1;
            return Index1D(j, k, type);
        }
    }

    std::cerr << "error: maptowav exit failure\n";
    exit(EXIT_FAILURE);
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_INDEXOPS_TCC
