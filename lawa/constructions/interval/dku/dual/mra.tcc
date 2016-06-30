/*
  This file is part of LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2011  Mario Rometsch, Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

#include <cassert>
#include <cmath>

#include <lawa/aux/aux.h>
#include <lawa/math/lawa_math.h>
#include <lawa/constructions/realline/dual/bspline.h>
#include <extensions/flens/lapack_flens.h>

namespace lawa {

template <typename T>
MRA<T,Dual,Interval,DKU>::MRA(FLENS_DEFAULT_INDEXTYPE _d, FLENS_DEFAULT_INDEXTYPE _d_, FLENS_DEFAULT_INDEXTYPE j)
    : d(_d), d_(_d_), mu(d&1),
      l1((mu-d)/2), l2((mu+d)/2),
      l1_(l1-d_+1), l2_(l2+d_-1),
      l_(l2_), q_(l_+mu-1),
      l(l_-(d_-d)), q(l+mu-1),
      min_j0(iceil<FLENS_DEFAULT_INDEXTYPE>(log2(l_+l2_-1)+1)),
      j0(std::max(j,min_j0)),
      _bc(2,0), _j(j0)
{
    assert(_j>=min_j0);
    assert((d+d_)%2==0);
    setLevel(_j);

    _alpha();
    _beta();
    
    _calcM0_();

    _integral0toInfPhiPhi_();
}

template <typename T>
MRA<T,Dual,Interval,DKU>::~MRA()
{
}

//--- cardinalities of whole, left, inner, right index sets. -------------------

template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Dual,Interval,DKU>::cardI_(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=min_j0);
    return d_ + pow2i<T>(j)-2*l_-mu+1 + d_ - (_bc(0)+_bc(1));
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Dual,Interval,DKU>::cardI_L(FLENS_DEFAULT_INDEXTYPE /*j*/) const
{
    return d_-_bc(0);
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Dual,Interval,DKU>::cardI_I(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=min_j0);
    return pow2i<T>(j)-2*l_-mu+1;
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Dual,Interval,DKU>::cardI_R(FLENS_DEFAULT_INDEXTYPE /*j*/) const
{
    return d_-_bc(1);
}

//--- ranges of whole, left, inner, right index sets. --------------------------

template <typename T>
flens::Range<FLENS_DEFAULT_INDEXTYPE>
MRA<T,Dual,Interval,DKU>::rangeI_(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=min_j0);
    return flens::Range<FLENS_DEFAULT_INDEXTYPE>(l_-d_+_bc(0), pow2i<T>(j)-l_-mu+d_-_bc(1));
}

template <typename T>
flens::Range<FLENS_DEFAULT_INDEXTYPE>
MRA<T,Dual,Interval,DKU>::rangeI_L(FLENS_DEFAULT_INDEXTYPE /*j*/) const
{
    return flens::Range<FLENS_DEFAULT_INDEXTYPE>(l_-d_+_bc(0),l_-1);
}

template <typename T>
flens::Range<FLENS_DEFAULT_INDEXTYPE>
MRA<T,Dual,Interval,DKU>::rangeI_I(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=min_j0);
    return flens::Range<FLENS_DEFAULT_INDEXTYPE>(l_, pow2i<T>(j)-l_-mu);
}

template <typename T>
flens::Range<FLENS_DEFAULT_INDEXTYPE>
MRA<T,Dual,Interval,DKU>::rangeI_R(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=min_j0);
    return flens::Range<FLENS_DEFAULT_INDEXTYPE>(pow2i<T>(j)-l_-mu+1, pow2i<T>(j)-l_-mu+d_-_bc(1));
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
MRA<T,Dual,Interval,DKU>::level() const
{
    return _j;
}

template <typename T>
void
MRA<T,Dual,Interval,DKU>::setLevel(FLENS_DEFAULT_INDEXTYPE j)
{
    if (j!=_j) {
        assert(j>=min_j0);
        _j = j;
    }
}

template <typename T>
template <BoundaryCondition BC>
void
MRA<T,Dual,Interval,DKU>::enforceBoundaryCondition()
{
    assert(BC==DirichletBC);

    _bc(0) = DirichletBC;
    _bc(1) = DirichletBC;
    
    _calcM0_();
}

template <typename T>
void
MRA<T,Dual,Interval,DKU>::_alpha()
{
    using flens::_;

    _Alpha.engine().resize(_(1-l2_, 2*l_+l1_-1), _(0,d_-1));

    _Alpha(_,0) = 1.;      // (5.1.1)

    // 5.1.3
     BSpline<T,Primal,R,CDF> phi(d);
     for (FLENS_DEFAULT_INDEXTYPE r=1; r<d_; ++r) {
         T tmp1 = 0.;
         for (FLENS_DEFAULT_INDEXTYPE k=l1; k<=l2; ++k) {
             T tmp2 = 0.;
             for (FLENS_DEFAULT_INDEXTYPE s=0; s<r; ++s) {
                 tmp2 += binomial(r,s)*powii(k,r-s)*_Alpha(0,s);
             }
             tmp1 += phi.a(k)*tmp2;
         }
        _Alpha(0,r) = tmp1/(pow2i<T>(r+1)-2);
    }

    // (5.1.2)
    for (FLENS_DEFAULT_INDEXTYPE r=1; r<d_; ++r) {
        // m<0
        for (FLENS_DEFAULT_INDEXTYPE m=_Alpha.firstRow(); m<0; ++m) {
            T tmp = 0.;
            for (FLENS_DEFAULT_INDEXTYPE i=0; i<=r; ++i) {
                tmp += binomial(r,i)*powii(m,i)*_Alpha(0,r-i);
            }
            _Alpha(m,r) = tmp;
        }
        // m>0
        for (FLENS_DEFAULT_INDEXTYPE m=1; m<=_Alpha.lastRow(); ++m) {
            T tmp = 0.;
            for (FLENS_DEFAULT_INDEXTYPE i=0; i<=r; ++i) {
                tmp += binomial(r,i)*powii(m,i)*_Alpha(0,r-i);
            }
            _Alpha(m,r) = tmp;
        }
    }
}

template <typename T>
void
MRA<T,Dual,Interval,DKU>::_beta()
{
    using flens::_;

    if (d==1) {
        _Beta.engine().resize(_(2*l_+l1_, 2*l_+l1_), _(0,0));
        return;
    }
    _Beta.engine().resize(_(2*l_+l1_, 2*l_+l2_-2), _(0,d_-1));

    // 3.2.31
    BSpline<T,Dual,R,CDF> phi_(d,d_);
    T factor = Const<T>::R_SQRT2;
    for (FLENS_DEFAULT_INDEXTYPE r=0; r<d_; ++r) {
        for (FLENS_DEFAULT_INDEXTYPE m=_Beta.firstRow(); m<=_Beta.lastRow(); ++m) {
            T tmp = 0.;
            for (FLENS_DEFAULT_INDEXTYPE q=iceil<FLENS_DEFAULT_INDEXTYPE>((m-l2_)/2.); q<l_; ++q) {
                tmp += _Alpha(q,r) * phi_.a_(m-2*q);
            }
            _Beta(m,r) = factor * tmp;
        }
    }
}

template <typename T>
void
MRA<T,Dual,Interval,DKU>::_calcM0_()
{
    using flens::_;

    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Mj0_(rangeI_(min_j0+1), rangeI_(min_j0));
    
    // left upper block
    FLENS_DEFAULT_INDEXTYPE upper = l2_+2*l_-2;
    FLENS_DEFAULT_INDEXTYPE low   = l_-d_;
    FLENS_DEFAULT_INDEXTYPE hi    = l_-1;
    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > ML_c(_(low,upper), _(low,hi)), MR_c;

    //     setup ML_c
    for (FLENS_DEFAULT_INDEXTYPE m=low; m<=hi; ++m) {
        ML_c(m,m) = Const<T>::R_SQRT2 * pow2i<T>(-m+low);
    }
    for (FLENS_DEFAULT_INDEXTYPE m=hi+1; m<=2*l_+l1_-1; ++m) {
        for (FLENS_DEFAULT_INDEXTYPE k=low; k<=hi; ++k) {
            ML_c(m,k) = Const<T>::R_SQRT2 * pow2i<T>(-k+low) * _Alpha(m,k-low);
        }
    }

    for (FLENS_DEFAULT_INDEXTYPE m=2*l_+l1_; m<=upper; ++m) {
        for (FLENS_DEFAULT_INDEXTYPE k=low; k<=hi; ++k) {
            ML_c(m,k) = _Beta(m,k-low);
        }
    }

    Mj0_(ML_c) = ML_c;
    arrow(ML_c,MR_c);
    MR_c.engine().changeIndexBase(Mj0_.lastRow()-MR_c.numRows()+1,
                                  Mj0_.lastCol()-MR_c.numCols()+1);
    Mj0_(MR_c) = MR_c;
    Mj0_.engine().changeIndexBase(1,1);

    // central band (vertical)
    BSpline<T,Dual,R,CDF> phi_(d,d_);
    for (FLENS_DEFAULT_INDEXTYPE c=l_, r=2*l_+l1_; c<=pow2i<T>(min_j0)-q_-1; ++c, r+=2) {
        Mj0_(_(r,r+phi_.a_.length()-1),c) = Const<T>::R_SQRT2 * phi_.a_; 
    }

    for (FLENS_DEFAULT_INDEXTYPE c=l_, r=2*l_+l1_; c<=pow2i<T>(min_j0)-q_-1; ++c, r+=2) {
        Mj0_(_(r,r+phi_.a_.length()-1),c) = Const<T>::R_SQRT2 * phi_.a_; 
    }
    
    M0_ = flens::RefinementMatrix<T,Interval,DKU>(d_,d_,Mj0_,min_j0);
}

template <typename T>
T
chi(T x)
{
    return (x>=0 && x<=1) ? 1 : 0;
}

template <typename T>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
MRA<T,Dual,Interval,DKU>::_integral0toInfPhiPhi_()
{
    using flens::_;

    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > I(_(-l2 +1, l-d+d_-1),
                                         _(-l2_+1, l_-1));
    BSpline<T,Primal,R,CDF> phi(d);
    BSpline<T,Dual,R,CDF> phi_(d,d_);

    return I;
}


} // namespace lawa

