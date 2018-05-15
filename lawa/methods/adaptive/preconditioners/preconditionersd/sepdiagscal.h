#ifndef LAWA_METHODS_ADAPTIVE_PRECONDITIONERS_PRECONDITIONERSD_SEPDIAGSCAL_H
#define LAWA_METHODS_ADAPTIVE_PRECONDITIONERS_PRECONDITIONERSD_SEPDIAGSCAL_H 1

#include <cstddef>
#include <vector>

#include <flens/flens.cxx>

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/mapwavind.h>

namespace lawa
{

template <typename Basis>
class Sepdiagscal
{
public:
    typedef std::size_t                            size_type;
    typedef FLENS_DEFAULT_INDEXTYPE                Int;
    typedef double                                 T;
    typedef flens::DenseVector<flens::Array<T> >   DenseVector;
    typedef flens::DenseVector<flens::Array<Int> > DenseIntVector;
    typedef flens::DiagMatrix<flens::Array<T> >    DiagMat;
    typedef std::vector<DiagMat>                   DiagMatVec;
    typedef Mapwavind<Index1D>                     Maptype;

private:
    const size_type     dim_;
    const Basis*        basis_;
          Maptype*      map_;
    const T             order_;
    DiagMatVec          Ds_;
    T                   eps_;
    T                   nu_;
    T                   h_;
    T                   iscale_;
    size_type           nplus_;
    size_type           n_;

public:
    Sepdiagscal()                   = delete;

    Sepdiagscal(const Sepdiagscal&) = default;

    Sepdiagscal(Sepdiagscal&&)      = default;

    Sepdiagscal(const size_type _dim, const Basis& _basis, Maptype& _map,
                const T _order = 1,
                const T _eps = 1., const T _nu = 1., const T _h = 0.5,
                const T _iscale = 2.,
                const size_type _nplus = 1, const size_type _n = 1);

    size_type
    dim() const;

    const Basis&
    basis() const;

    Maptype&
    map();

    T
    order() const;

    T
    eps() const;

    T
    nu() const;

    T
    h() const;

    T
    iscale() const;

    size_type
    nplus() const;

    size_type
    n() const;

    void
    set_eps(const T _eps);

    void
    set_nu(const T _nu);

    void
    set_h(const T _h);

    void
    comp_h();

    void
    set_iscale(const T _iscale);

    void
    set_nplus(const size_type _nplus);

    void
    comp_nplus();

    void
    set_n(const size_type _n);

    void
    comp_n();

    void
    assemble(const std::vector<IndexSet<Index1D> >& cols);

    const DiagMat&
    operator()(const Int k, const size_type j) const;

    DiagMat&
    operator()(const Int k, const size_type j);
};

} // namespace lawa

#include <lawa/methods/adaptive/preconditioners/preconditionersd/sepdiagscal.tcc>

#endif // LAWA_METHODS_ADAPTIVE_PRECONDITIONERS_PRECONDITIONERSD_SEPDIAGSCAL_H
