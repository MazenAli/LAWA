#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_MAPWAVIND_TCC
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_MAPWAVIND_TCC 1

#include <cassert>
#include <stdlib.h>
#include <limits>
#include <unordered_map>

namespace lawa
{

template <typename Index>
Mapwavind<Index>::Mapwavind(const size_type d):
    activex(d)
{
    for (size_type i=0; i<d; ++i) {
        activex[i].left.rehash(DEFAULT_HASH_SIZE);
        activex[i].right.rehash(DEFAULT_HASH_SIZE);
    }
}


template <typename Index>
typename Mapwavind<Index>::size_type
Mapwavind<Index>::dim() const
{
    return activex.size();
}


template <typename Index>
const typename Mapwavind<Index>::IndexSetVec&
Mapwavind<Index>::get_active() const
{
    return activex;
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

    activex[d-1].left.rehash(n);
    activex[d-1].right.rehash(n);
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
typename Mapwavind<Index>::Uint
Mapwavind<Index>::operator()(const Index& index, const size_type d)
{
    assert(d>=1 && d<=dim());

    /* Retrive index or update activex */
    auto loc  = activex[d-1].left.find(index);
    Uint val  = 0;
    if (loc==activex[d-1].left.end()) {
        val = activex[d-1].left.size()+1;
        std::pair<Index, Uint> newind(index, val);
        activex[d-1].left.insert(newind);
    } else {
        val = loc->second;
    }

    return val;
}


template <typename Index>
Index
Mapwavind<Index>::operator()(const Uint i,
                             const size_type d) const
{
    assert(d>=1 && d<=dim());

    /* Retrive index or update activex */
    auto loc  = activex[d-1].right.find(i);
    if (loc!=activex[d-1].right.end()) {
        return loc->second;
    }

    std::cerr <<
    "error mapwavind::operator(): wavelet index has not been mapped before\n";

    return Index();
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_MAPWAVIND_TCC
