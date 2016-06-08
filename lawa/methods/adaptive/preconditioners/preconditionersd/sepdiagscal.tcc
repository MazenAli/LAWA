#ifndef LAWA_METHODS_ADAPTIVE_PRECONDITIONERS_PRECONDITIONERSD_SEPDIAGSCAL_TCC
#define LAWA_METHODS_ADAPTIVE_PRECONDITIONERS_PRECONDITIONERSD_SEPDIAGSCAL_TCC 1

#define _USE_MATH_DEFINES
#ifndef MAX
    #define MAX(x, y) (x>y) ? x : y
#endif
#ifndef MIN
    #define MIN(x, y) (x<y) ? x : y
#endif

#include <cassert>
#include <cmath>

namespace lawa
{

template <typename Basis>
Sepdiagscal<Basis>::Sepdiagscal(const size_type _dim, const Basis& _basis,
                                const T _order, const T _eps, const T _nu,
                                const T _h, const T _iscale,
                                const size_type _nplus, const size_type _n):
    dim_(_dim),
    basis_(&_basis),
    order_(_order),
    eps_(_eps),
    nu_(_nu),
    h_(_h),
    iscale_(_iscale),
    nplus_(_nplus),
    n_(_n)
{
    assert(dim()>0);
    assert(eps()>0);
    assert(nu()>0);
    assert(h()>0);
    assert(iscale()>0);
    assert(nplus()>0);
    assert(n()>0);
}


template <typename Basis>
typename Sepdiagscal<Basis>::size_type
Sepdiagscal<Basis>::dim() const
{
    return dim_;
}


template <typename Basis>
const Basis&
Sepdiagscal<Basis>::basis() const
{
    assert(basis_);
    return *basis_;
}


template <typename Basis>
typename Sepdiagscal<Basis>::T
Sepdiagscal<Basis>::order() const
{
    return order_;
}


template <typename Basis>
typename Sepdiagscal<Basis>::T
Sepdiagscal<Basis>::eps() const
{
    return eps_;
}


template <typename Basis>
typename Sepdiagscal<Basis>::T
Sepdiagscal<Basis>::nu() const
{
    return nu_;
}


template <typename Basis>
typename Sepdiagscal<Basis>::T
Sepdiagscal<Basis>::h() const
{
    return h_;
}


template <typename Basis>
typename Sepdiagscal<Basis>::T
Sepdiagscal<Basis>::iscale() const
{
    return iscale_;
}


template <typename Basis>
typename Sepdiagscal<Basis>::size_type
Sepdiagscal<Basis>::nplus() const
{
    return nplus_;
}


template <typename Basis>
typename Sepdiagscal<Basis>::size_type
Sepdiagscal<Basis>::n() const
{
    return n_;
}


template <typename Basis>
void
Sepdiagscal<Basis>::set_eps(const T _eps)
{
    assert(_eps>0);
    eps_ = _eps;
}


template <typename Basis>
void
Sepdiagscal<Basis>::set_nu(const T _nu)
{
    assert(_nu>0);
    nu_ = _nu;
}


template <typename Basis>
void
Sepdiagscal<Basis>::set_h(const T _h)
{
    assert(_h>0);
    h_ = _h;
}


template <typename Basis>
void
Sepdiagscal<Basis>::comp_h()
{
    T h_ = (M_PI*M_PI)/(5.*(std::abs(std::log(eps()/2.))+4.));
    set_h(h_);
}


template <typename Basis>
void
Sepdiagscal<Basis>::set_iscale(const T _iscale)
{
    assert(_iscale>0.);
    iscale_ = _iscale;
}


template <typename Basis>
void
Sepdiagscal<Basis>::set_nplus(const size_type _nplus)
{
    assert(_nplus>0);
    nplus_ = _nplus;
}


template <typename Basis>
void
Sepdiagscal<Basis>::comp_nplus()
{
    T max            = MAX(4.*std::pow(M_PI, -0.5),
                           std::sqrt(std::abs(std::log(eps()/2.))));
    size_type _nplus = std::ceil((1./h())*max);
    set_nplus(_nplus);
}


template <typename Basis>
void
Sepdiagscal<Basis>::set_n(const size_type _n)
{
    assert(_n>0);
    n_ = _n;
}


template <typename Basis>
void
Sepdiagscal<Basis>::comp_n()
{
    T min        = MIN(eps()/2., nu());
    size_type _n = std::ceil((1./h())*(std::log(2.*std::pow(M_PI,-0.5))+
                   std::abs(std::log(min))
                   +0.5*std::log(iscale())));
    set_n(_n);
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_PRECONDITIONERS_PRECONDITIONERSD_SEPDIAGSCAL_TCC 1
