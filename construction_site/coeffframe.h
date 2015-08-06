#include <cstddef>
#include <vector>
#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>

#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_COEFFFRAME_H
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_COEFFFRAME_H 1

namespace lawa
{


template <SortingCriterion S, typename T, typename _Index>
class CoeffFrame
{
    private:
        const std::size_t                           numCols_;
        std::vector<Coefficients<S, T, _Index> >    U;
    public:
        // types
        typedef typename Coefficients<S, T, _Index>::const_iterator  const_iterator;
        typedef typename Coefficients<S, T, _Index>::iterator        iterator;
        typedef _Index                                               IndexType;

        // construct/copy/destruct
        CoeffFrame(const std::size_t _numCols_);

        CoeffFrame(const CoeffFrame<S, T, _Index>& copy);

        // size
        std::size_t
        numCols() const;

        // observers
        const std::vector<Coefficients<S, T, _Index> >&
        getFrame() const;

        const Coefficients<S, T, _Index>&
        operator[] (const std::size_t col_num) const;

        T
        operator() (const _Index& lambda, const std::size_t col_num) const;

        // modifiers
        Coefficients<S, T, _Index>&
        operator[] (const std::size_t col_num);

        T&
        operator() (const _Index& lambda, const std::size_t col_num);
};


} // namespace lawa

#include <construction_site/coeffframe.tcc>

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_COEFFFRAME_H
