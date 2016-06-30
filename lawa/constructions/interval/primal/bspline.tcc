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
#include <algorithm>
#include <limits>

namespace lawa {

template <typename T>
    T
    _evaluateUnitBSpline(FLENS_DEFAULT_INDEXTYPE d, T x, FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k, unsigned short deriv);

//------------------------------------------------------------------------------

template <typename T, Construction Cons>
BSpline<T,Primal,Interval,Cons>::BSpline(const MRA<T,Primal,Interval,Cons> &_mra)
    : mra(_mra)
{
}

template <typename T, Construction Cons>
T
BSpline<T,Primal,Interval,Cons>::operator()(T x, FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k, unsigned short deriv) const
{
    assert(j>=mra.j0);
    assert(k>=mra.rangeI(j).firstIndex());
    assert(k<=mra.rangeI(j).lastIndex());
    return _evaluateUnitBSpline(mra.d, x, j, k, deriv);
}

template <typename T, Construction Cons>
Support<T>
BSpline<T,Primal,Interval,Cons>::support(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    assert(j>=mra.j0);
    assert(k>=mra.rangeI(j).firstIndex());
    assert(k<=mra.rangeI(j).lastIndex());
    return pow2i<T>(-j) * Support<T>(std::max((FLENS_DEFAULT_INDEXTYPE) 0,k-mra.d),
                                     std::min(k,pow2i<FLENS_DEFAULT_INDEXTYPE>(j)));
}

template <typename T, Construction Cons>
flens::DenseVector<flens::Array<T> >
BSpline<T,Primal,Interval,Cons>::singularSupport(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    const FLENS_DEFAULT_INDEXTYPE tics = (k<mra.d) ? k+1 : (k>pow2i<T>(j)) ? pow2i<T>(j)+mra.d-1-k+2 : mra.d+1;
    return linspace(pow2i<T>(-j) * std::max((FLENS_DEFAULT_INDEXTYPE) 0,k-mra.d),
                    pow2i<T>(-j) * std::min(T(k),pow2i<T>(j)),
                    tics);
}

template <typename T, Construction Cons>
T
BSpline<T,Primal,Interval,Cons>::tic(FLENS_DEFAULT_INDEXTYPE j) const
{
    return pow2i<T>(-j);
}

template <typename T, Construction Cons>
FLENS_DEFAULT_INDEXTYPE
BSpline<T,Primal,Interval,Cons>::getRefinementLevel(FLENS_DEFAULT_INDEXTYPE j) const
{
    return j + 1;
}

template <typename T, Construction Cons>
flens::DenseVector<flens::Array<long double> > *
BSpline<T,Primal,Interval,Cons>::getRefinement(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k, FLENS_DEFAULT_INDEXTYPE &refinement_j, FLENS_DEFAULT_INDEXTYPE &refinement_k_first,
												FLENS_DEFAULT_INDEXTYPE &split, FLENS_DEFAULT_INDEXTYPE &refinement_k_restart) const
{
	// No split necessary, so set default values
	refinement_k_restart = 1;

    refinement_j = j + 1;
    // left boundary
    if (k<mra.rangeII(j).firstIndex()) {
        FLENS_DEFAULT_INDEXTYPE type  = k % mra.cardIL(j);
        refinement_k_first = mra._leftOffsets[type];
        split = mra._leftRefCoeffs[type].length()+1;
        return &(mra._leftRefCoeffs[type]);
    }
    // inner part
    if (k<=mra.rangeII(j).lastIndex()) {
        refinement_k_first = 2*k+mra._innerOffsets[0];
        split = mra._innerRefCoeffs[0].length()+1;
        return &(mra._innerRefCoeffs[0]);
    }
    // right part
    //FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)(k - (mra.cardI(j)-1 - mra._numRightParts + 1));
    FLENS_DEFAULT_INDEXTYPE type  = (mra.rangeI(j).lastIndex()-k);
    FLENS_DEFAULT_INDEXTYPE shift = pow2i<FLENS_DEFAULT_INDEXTYPE>(j)-1;
    refinement_k_first =2*shift+mra._rightOffsets[type];
    split = mra._rightRefCoeffs[type].length()+1;
    return &(mra._rightRefCoeffs[type]);
}

template <typename T, Construction Cons>
T
BSpline<T,Primal,Interval,Cons>::getL2Norm(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    if (k<mra.rangeII(j).firstIndex()) {
        FLENS_DEFAULT_INDEXTYPE type  = k % mra.cardIL(j);
        return mra._leftL2Norms[type];
    }
    // inner part
    if (k<=mra.rangeII(j).lastIndex()) {
        return mra._innerL2Norms[0];
    }
    // right part
    FLENS_DEFAULT_INDEXTYPE type  = (mra.rangeI(j).lastIndex()-k);
    return mra._rightL2Norms[type];
}

template <typename T, Construction Cons>
T
BSpline<T,Primal,Interval,Cons>::getH1SemiNorm(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    long double pow2ij = (long double)((FLENS_DEFAULT_INDEXTYPE) 1 << j);
    if (k<mra.rangeII(j).firstIndex()) {
        FLENS_DEFAULT_INDEXTYPE type  = k % mra.cardIL(j);
        return pow2ij*mra._leftH1SemiNorms[type];
    }
    // inner part
    if (k<=mra.rangeII(j).lastIndex()) {
        return pow2ij*mra._innerH1SemiNorms[0];
    }
    // right part
    FLENS_DEFAULT_INDEXTYPE type  = (mra.rangeI(j).lastIndex()-k);
    return pow2ij*mra._rightH1SemiNorms[type];
}

//--- evaluate B-spline --------------------------------------------------------

template <typename T>
T
_evaluateUnitBSpline(FLENS_DEFAULT_INDEXTYPE d, T x, FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k, unsigned short deriv)
{
    assert(x>=0.0);
    assert(x<=1.0);
    
    if (deriv>=d) {
        return 0;
    }
    
    if (deriv==0) {
        // "if" needed for calculation of derivatives
        //  (otherwise assertion would be correct).
        if ((k>pow2i<T>(j)+d-1) || (k<1)) {
            return 0;
        }

        T twoj = pow2i<T>(j); 
        x *= twoj;
        FLENS_DEFAULT_INDEXTYPE pos = ifloor(x) - (x==twoj);
        // we are not inside the support.
        if ((pos<k-d) || (pos>k-1)) {
            return 0;
        }

        flens::Array<T> values(d,0);
        // initialize correct 'slot'.
        if ((pos<=k) && (pos>=k-d)) {
            values(pos-(k-d)) = 1.;
        }

        // utilizing left multiplicities.
        if (k<d) {
            for (FLENS_DEFAULT_INDEXTYPE m=2; m<=d; ++m) {
                for (FLENS_DEFAULT_INDEXTYPE i=0; i<=d-m; ++i) {
                    if (m+i+k-1>d) {
                        values(i) = (x-std::max((FLENS_DEFAULT_INDEXTYPE) 0,i+k-d))*values(i) / (m+i+k-1-d-std::max(k+i-d,(FLENS_DEFAULT_INDEXTYPE) 0));                    
                    }
                    if (m+i+k>d) {
                        values(i) += (m+i+k-d-x)*values(i+1) / (m+i+k-d-std::max(k+i+1-d,(FLENS_DEFAULT_INDEXTYPE) 0));                    
                    }
                }
            }
            return pow2ih<T>(j)*values(0);
        }

        // utilizing right multiplicities.
        FLENS_DEFAULT_INDEXTYPE t = twoj;
        if (k>t) {
            for (FLENS_DEFAULT_INDEXTYPE m=2; m<=d; ++m) {
                for (FLENS_DEFAULT_INDEXTYPE i=0; i<=d-m; ++i) {
                    if (k+i-d<t) {
                        values(i) = (x-(i+k-d))*values(i) / (std::min(m+i+k-1-d,t)-(k+i-d));                    
                    }
                    if (k+i+1-d<t) {
                        values(i) += (std::min(t,m+i+k-d)-x)*values(i+1) / (std::min(t,m+i+k-d)-(k+i+1-d));                    
                    }
                }
            }
            return pow2ih<T>(j)*values(0);
        }

        // 'inner' B-Splines.
        t = k-d;
        for (FLENS_DEFAULT_INDEXTYPE m=2; m<=d; ++m) {
            for (FLENS_DEFAULT_INDEXTYPE i=0; i<=d-m; ++i) {
                values(i) =  ((x-(t+i))*values(i) + ((t+m+i)-x)*values(i+1))/(m-1);
            }
        }
        return pow2ih<T>(j)*values(0);
    } else {
        assert(k>=1);
        assert(k<=pow2i<T>(j)+d-1);

        T value = 0.;

        T twomj = pow2i<T>(-j);
        T a=twomj*(k-1), b=twomj*(k), c=twomj*(k-d), e=twomj*(k-d+1);
        if (k<=d) {
            c = 0;
        }
        if (k+d>=pow2i<T>(j)+d) {
            b = 1.;
        }
        if (k+d-1<=d) {
            a = 0.;
        }
        if (k+d-1>=pow2i<T>(j)+d) {
            a = 1.;
        }
        if (k+1<=d) {
            e = 0.;
        }

        // remark: the index k is shifted by -1 since for d-1 the implicit knot
        //         vector is implicitely shifted to the left by one.
        if (a!=c) {
            value  = _evaluateUnitBSpline(d-1,x,j,k-1,deriv-1)   / (a-c);                    
        }
        if (b!=e) {
            value -= _evaluateUnitBSpline(d-1,x,j,k,deriv-1) / (b-e);                    
        }

        return (d-1)*value;
    }
}

} // namespace lawa

