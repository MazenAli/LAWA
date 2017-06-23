#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_MAPWAVIND_H
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_MAPWAVIND_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <flens/flens.cxx>
#include <vector>

#include <boost/bimap.hpp>
#include <boost/bimap/unordered_set_of.hpp>

#define DEFAULT_HASH_SIZE 1000

namespace lawa
{

template <typename Index>
class Mapwavind
{
public:
    typedef unsigned FLENS_DEFAULT_INDEXTYPE              Uint;
    typedef          boost::bimaps::unordered_set_of<
                     Index,
                     index_hashfunction<Index>,
                     index_eqfunction<Index>>             LeftType;
    typedef          boost::bimaps::unordered_set_of<
                     Uint>                                RightType;
    typedef typename boost::bimap<LeftType, RightType>    Maptype;
    typedef typename std::vector<Maptype>                 IndexSetVec;
    typedef typename IndexSetVec::size_type               size_type;

private:
    IndexSetVec  activex;

public:
    Mapwavind()                 = default;
    Mapwavind(const Mapwavind&) = default;
    Mapwavind(Mapwavind&&)      = default;

    Mapwavind(const size_type d);

    size_type
    dim() const;

    const IndexSetVec&
    get_active() const;

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

    Uint
    operator()(const Index& index, const size_type d);

    Index
    operator()(const Uint i,
               const size_type d) const;
};

} // namespace lawa

#include <lawa/methods/adaptive/datastructures/mapwavind.tcc>

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_MAPWAVIND_H
