#ifndef LAWA_METHODS_ADAPTIVE_OPERATORS_OPERATORSD_SEPOP_H
#define LAWA_METHODS_ADAPTIVE_OPERATORS_OPERATORSD_SEPOP_H 1

#include <vector>
#include <cstddef>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/htcoefficients.h>
#include <lawa/settings/enum.h>

namespace lawa
{

template <typename Optype>
class Sepop
{
public:
    typedef typename std::vector<Optype*>               Opvec;
    typedef typename std::vector<Optype*>::size_type    size_type;

private:
    Opvec               ops_;
    const size_type     rank_;
    const size_type     dim_;
    IndexSet<Index1D>   indexset_;
    SepopType           type_;

public:
    Sepop()             = delete;

    Sepop(const Sepop&) = default;

    Sepop(Sepop&&)      = default;

    Sepop(const Opvec& _ops, const size_type _rank, const size_type _dim);

    Sepop(const Optype& op, const size_type _rank, const size_type _dim);

    size_type
    rank() const;

    size_type
    dim() const;

    SepopType
    type() const;

    const IndexSet<Index1D>&
    getIndexset() const;

    void
    setIndexset(const IndexSet<Index1D>& _indexset);

    const Opvec&
    ops() const;

    Opvec&
    ops();

    const Optype&
    ops(const size_type i, const size_type j) const;

    Optype&
    ops(const size_type i, const size_type j);

    const Optype&
    operator()(const size_type i, const size_type j) const;

    Optype&
    operator()(const size_type i, const size_type j);
/*
    template <typename T, typename Basis>
    HTCoefficients<T, Basis>
    eval(const HTCoefficients<T, Basis>& u,
         const IndexSet<Index1D>& indexset,
         const std::size_t hashtablelength=193);
*/
};

} // namespace lawa

#include <lawa/methods/adaptive/operators/operatorsd/sepop.tcc>

#endif // LAWA_METHODS_ADAPTIVE_OPERATORS_OPERATORSD_SEPOP_H
