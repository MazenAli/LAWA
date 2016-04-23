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
    typedef typename std::vector<IndexSet<Index1D> >    IndexSetVec;
    typedef typename std::vector<Optype*>::size_type    size_type;

private:
    Opvec               ops_;
    const size_type     rank_;
    const size_type     dim_;
    IndexSetVec         indrows;
    IndexSetVec         indcols;
    SepopType           type_;

public:
    Sepop()             = delete;

    Sepop(const Sepop&) = default;

    Sepop(Sepop&&)      = default;

    Sepop(const Opvec& _ops, const size_type _rank, const size_type _dim);

    Sepop(Optype& op, const size_type _rank, const size_type _dim);

    size_type
    rank() const;

    size_type
    dim() const;

    SepopType
    type() const;

    const IndexSetVec&
    getrows() const;

    const IndexSetVec&
    getcols() const;

    const IndexSet<Index1D>&
    getrows(const size_type j) const;

    const IndexSet<Index1D>&
    getcols(const size_type j) const;

    void
    setrows(const IndexSetVec& _indrows);

    void
    setcols(const IndexSetVec& _indcols);

    void
    setrows(const IndexSet<Index1D>& _indrows, const size_type j);

    void
    setcols(const IndexSet<Index1D>& _indcols, const size_type j);


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
};

} // namespace lawa

#include <lawa/methods/adaptive/operators/operatorsd/sepop.tcc>

#endif // LAWA_METHODS_ADAPTIVE_OPERATORS_OPERATORSD_SEPOP_H
