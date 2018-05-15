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
#include <lawa/methods/adaptive/algorithms/coeffops.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>

namespace lawa
{

template <typename Basis>
Sepdiagscal<Basis>::Sepdiagscal(const size_type _dim, const Basis& _basis,
                                      Maptype& _map,
                                const T _order, const T _eps, const T _nu,
                                const T _h, const T _iscale,
                                const size_type _nplus, const size_type _n):
    dim_(_dim),
    basis_(&_basis),
    map_(&_map),
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
typename Sepdiagscal<Basis>::Maptype&
Sepdiagscal<Basis>::map()
{
    assert(map_);
    return *map_;
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
    T max            = std::max(4./std::sqrt(M_PI),
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


template <typename Basis>
void
Sepdiagscal<Basis>::assemble(const std::vector<IndexSet<Index1D> >& cols)
{
    assert(cols.size() == dim());

    // n and n+
    T iscale = compIndexscale(basis(), cols, order());
    set_iscale(iscale);
    comp_n();
    size_type rank = n() + nplus() + 1;

    // Set up storage
    Ds_.resize(dim()*rank);

    // w_min
    T omega2  = compOmegamin2(basis(), dim(), order());
    T factor1 = h()*(1./std::sqrt(omega2));

    DenseIntVector  sizes(dim());
    for (size_type j=1; j<=dim(); ++j) {
        sizes(j) = maxintindhash(cols[j-1], j, map());
    }

    for (Int i=-1*n(); i<=(signed) nplus(); ++i) {
        T factor2 = 2.*std::pow(M_PI, -0.5)*
                    (1./(1+std::exp(-1.*h()*(T) i)));
        factor2   = std::pow(factor1*factor2, 1./(T) dim());
        T alpha   = std::pow(std::log(1.+std::exp((T) i*h())), 2.);

        for (size_type j=1; j<=dim(); ++j) {
            DenseVector Dij(sizes(j));

            for (const auto& lambda : cols[j-1]) {
                Int level = lambda.j;
                if (lambda.xtype==XWavelet) ++level;
                T weight  = std::pow(2., 2.*order()*level)/omega2;

                auto k    = map()(lambda, j);
                Dij(k) = factor2*std::exp(-alpha*weight);
            }
            operator()(i, j).resize(sizes(j));
            operator()(i, j).diag() = Dij;
        }
    }
}


template <typename Basis>
const typename Sepdiagscal<Basis>::DiagMat&
Sepdiagscal<Basis>::operator()(const Int       k,
                               const size_type j) const
{
    assert(k>=-1*(signed) n() && k<=(signed) nplus());
    const size_type i    = k+(signed) n();
    const size_type rank = nplus()+n()+1;
    assert(j>=1 && j <= dim());
    assert((j-1)*rank+i <= Ds_.size());

    return Ds_[(j-1)*rank+i];
}


template <typename Basis>
typename Sepdiagscal<Basis>::DiagMat&
Sepdiagscal<Basis>::operator()(const Int       k,
                               const size_type j)
{
    assert(k>=-1*(signed) n() && k<=(signed) nplus());
    const size_type i    = k+(signed) n();
    const size_type rank = nplus()+n()+1;
    assert(j>=1 && j <= dim());
    assert((j-1)*rank+i <= Ds_.size());

    return Ds_[(j-1)*rank+i];
}

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_PRECONDITIONERS_PRECONDITIONERSD_SEPDIAGSCAL_TCC 1
