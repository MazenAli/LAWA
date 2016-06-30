/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

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
#include <lawa/flensforlawa.h>

namespace flens {

template <typename T>
template <lawa::FunctionSide Side>
RefinementMatrix<T,lawa::R,lawa::CDF>::RefinementMatrix(const lawa::BSpline<T,Side,lawa::R,lawa::CDF> &spline)
      : band(lawa::Const<T>::R_SQRT2 * spline.mask())
{
}

template <typename T>
template <lawa::FunctionSide Side>
RefinementMatrix<T,lawa::R,lawa::CDF>::RefinementMatrix(const lawa::Wavelet<T,Side,lawa::R,lawa::CDF> &wavelet)
      : band(lawa::Const<T>::R_SQRT2 * wavelet.mask())
{
}


//------------------------------------------------------------------------------

/*
template <typename X, typename Y>
void
mv(cxxblas::Transpose transA, typename X::ElementType alpha,
   const RefinementMatrix<typename X::ElementType,lawa::R,lawa::CDF> &A,
   const flens::DenseVector<X> &x, typename X::ElementType beta, flens::DenseVector<Y> &y)
{
    assert(0);
    typedef typename X::ElementType T;

    assert(alpha==1.);
    assert(x.engine().stride()==1);

    const flens::DenseVector<flens::Array<T> > &a = A.band;
    FLENS_DEFAULT_INDEXTYPE lx = x.length();
    FLENS_DEFAULT_INDEXTYPE la = a.length();
    
    if (transA==cxxblas::NoTrans) {
        if (beta==0) {
            y.engine().resize(2*lx,x.firstIndex()-a.firstIndex()) || y.engine().fill();
        } else {
            assert(y.length()==2*lx);
            y.engine().changeIndexBase(x.firstIndex()-a.firstIndex());
        }
        const T *xp = x.engine().data();
        for (FLENS_DEFAULT_INDEXTYPE k=0; k<lx; ++k, ++xp) {
            FLENS_DEFAULT_INDEXTYPE mMin = a.firstIndex() + 2*k;
            const T *abegin = a.engine().data();
                  T *ybegin = y.engine().data();
            cxxblas::axpy(la, *xp, abegin, 1, ybegin+mMin, 1);
        }
    } else { // (transA==cxxblas::Trans)
        if (beta==0) {
            y.engine().resize(lx/2,x.firstIndex()+a.firstIndex()) || y.engine().fill();
        } else {
            assert(y.length()==lx/2);
            y.engine().changeIndexBase(x.firstIndex()+a.firstIndex());
        }
        T *iter = y.engine().data();
        for (FLENS_DEFAULT_INDEXTYPE m=0; m<lx/2; ++m, ++iter) {
            FLENS_DEFAULT_INDEXTYPE kMin = 2*m;
            const T *abegin = a.engine().data();
            T dotValue;
            cxxblas::dot(la, abegin, 1, x.engine().data()+kMin, 1, dotValue);
            *iter += dotValue;
        }
    }
}
*/

} // namespace flens

