#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_MAPWAVIND_H
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_MAPWAVIND_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <flens/flens.cxx>
#include <vector>

#define DEFAULT_HASH_SIZE 1000

namespace lawa
{

template <typename Index>
class Mapwavind
{
public:
    typedef std::vector<IndexSet<Index> >   IndexSetVec;
    typedef typename IndexSetVec::size_type size_type;

private:
    IndexSetVec  activex;

public:
    Mapwavind()                 = default;
    Mapwavind(const Mapwavind&) = default;
    Mapwavind(Mapwavind&&)      = default;

    Mapwavind(const size_type d);

    size_type
    dim() const;

    IndexSetVec&
    get_active();

    void
    rehash(const size_type n, const size_type d);

    void
    rehash(const size_type n);

    size_type
    size(const size_type d);

    size_type
    buckets(const size_type d);

    double
    load(const size_type d);

    unsigned FLENS_DEFAULT_INDEXTYPE
    operator()(const Index& index, const size_type d);

    Index
    operator()(const unsigned FLENS_DEFAULT_INDEXTYPE i,
               const size_type d) const;
};

} // namespace lawa

#include <lawa/methods/adaptive/datastructures/mapwavind.tcc>

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_MAPWAVIND_H
