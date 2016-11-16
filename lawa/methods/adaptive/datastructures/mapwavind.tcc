#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_MAPWAVIND_TCC
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_MAPWAVIND_TCC 1

#include <cassert>
#include <stdlib.h>
#include <limits>

namespace lawa
{

template <typename Index>
Mapwavind<Index>::Mapwavind(const size_type d):
    activex(d)
{
    for (size_type i=0; i<d; ++i) {
        activex[i].max_load_factor(std::numeric_limits<double>::infinity());
        activex[i].rehash(DEFAULT_HASH_SIZE);
    }
}


template <typename Index>
typename Mapwavind<Index>::size_type
Mapwavind<Index>::dim() const
{
    return activex.size();
}


template <typename Index>
typename Mapwavind<Index>::IndexSetVec&
Mapwavind<Index>::get_active()
{
    return activex;
}


template <typename Index>
void
Mapwavind<Index>::rehash(const size_type n, const size_type d)
{
    assert(d>=1 && d<=dim());

    activex[d-1].rehash(n);
}


template <typename Index>
void
Mapwavind<Index>::rehash(const size_type n)
{
    for (size_type i=1; i<=dim(); ++i) {
        rehash(n, i);
    }
}


template <typename Index>
typename Mapwavind<Index>::size_type
Mapwavind<Index>::size(const size_type d)
{
    assert(d>=1 && d<=dim());

    return activex[d-1].size();
}


template <typename Index>
typename Mapwavind<Index>::size_type
Mapwavind<Index>::buckets(const size_type d)
{
    assert(d>=1 && d<=dim());

    return activex[d-1].bucket_count();
}


template <typename Index>
double
Mapwavind<Index>::load(const size_type d)
{
    assert(d>=1 && d<=dim());

    return activex[d-1].load_factor();
}


template <typename Index>
unsigned FLENS_DEFAULT_INDEXTYPE
Mapwavind<Index>::operator()(const Index& index, const size_type d)
{
    assert(d>=1 && d<=dim());

    unsigned FLENS_DEFAULT_INDEXTYPE    ret = 0, count = 1;
    size_type                           size, size_local;
    index_eqfunction<Index>             eq;

    /* Update activex */
    if (activex[d-1].find(index)==activex[d-1].cend()) {
        activex[d-1].insert(index);
    }

    ret  = activex[d-1].bucket(index);
    size = activex[d-1].bucket_count();

    /* Count collisions */
    for (auto it =activex[d-1].cbegin(ret);
              it!=activex[d-1].cend(ret); ++it, ++count) {
        if (eq(*it, index)) break;
    }
    size_local = activex[d-1].bucket_size(activex[d-1].bucket(index));

    return ret+1+size*(size_local-count);
}


template <typename Index>
Index
Mapwavind<Index>::operator()(const unsigned FLENS_DEFAULT_INDEXTYPE i,
                             const size_type d) const
{
    unsigned FLENS_DEFAULT_INDEXTYPE n1, n2, count = 1;
    size_type size = activex[d].bucket_count(), size_local;

    n2 = i/size;
    n1 = i%size-1;

    if (n1>=size || n2>=activex[d-1].bucket_size(n1)) {
        std::cerr <<
        "error mapwavind::operator(): wavelet index has not been mapped before\n";
        exit(1);
    }

    size_local = activex[d-1].bucket_size(n1);
    for (auto it = activex[d-1].cbegin(n1);
              it!= activex[d-1].cend(n1); ++it, ++count) {
        if ((size_local-count)==n2) return *it;
    }

    std::cerr << "error maptowavind::operator(): unexpected error\n";
    exit(1);
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_MAPWAVIND_TCC
