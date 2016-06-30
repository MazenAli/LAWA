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
#include <lawa/flensforlawa.h>

#include <lawa/math/lawa_math.h>
#include <lawa/constructions/realline/primal/bspline.h>
#include <extensions/extensions.h>

namespace lawa {

//template <typename T>
//BSpline<T,Primal,Periodic,CDF>::BSpline(FLENS_DEFAULT_INDEXTYPE _d)
//    : d(_d), mu(d&1), phiR(_d)
//{
//    assert(_d>0);
//}

template <typename T>
BSpline<T,Primal,Periodic,CDF>::BSpline(const MRA<T,Primal,Periodic,CDF> &_mra)
    : d(_mra.d), mu(d&1), phiR(d), mra(_mra)
{
    assert(d>0);
}

template <typename T>
BSpline<T,Primal,Periodic,CDF>::~BSpline()
{
}

template <typename T>
T
BSpline<T,Primal,Periodic,CDF>::operator()(T x, FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k, unsigned short deriv) const
{
    // maximal support: [0,1]
    if((x < 0.) || (x > 1.)){
        return 0.;
    }
    
    // add contributions of original spline on R
    // = 'wrapping' around [0,1]
    T val = 0;
    for(FLENS_DEFAULT_INDEXTYPE l = ifloor(phiR.support(j,k).l1); l < iceil<FLENS_DEFAULT_INDEXTYPE>(phiR.support(j,k).l2); ++l){
        val += phiR(l+x, j, k, deriv);
    }
    return val;
    
}

template <typename T>
PeriodicSupport<T>
BSpline<T,Primal,Periodic,CDF>::support(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    Support<T> suppR = phiR.support(j,k);
    if(suppR.length() >= 1){
        return PeriodicSupport<T>(0,1);
    }
    if(suppR.l1 < 0){
        return PeriodicSupport<T>(0,1,suppR.l2, suppR.l1 + 1);
    }
    if(suppR.l2 > 1){
        return PeriodicSupport<T>(0,1,suppR.l2 - 1, suppR.l1);
    }
    return PeriodicSupport<T>(suppR.l1, suppR.l2);
}


template <typename T>
flens::DenseVector<flens::Array<T> >
BSpline<T,Primal,Periodic,CDF>::singularSupport(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{   
    if((phiR.support(j,k).l1 >= 0) && (phiR.support(j,k).l2 <= 1)){
         return linspace(support(j,k).l1, support(j,k).l2, d+1);
    }
    
    std::list<T> temp;
    flens::DenseVector<flens::Array<T> > singSuppR = linspace(phiR.support(j,k).l1, phiR.support(j,k).l2, d+1);
    temp.push_back(0.);
    temp.push_back(1.);
    for(FLENS_DEFAULT_INDEXTYPE i = singSuppR.firstIndex(); i <= singSuppR.lastIndex(); ++i){
        temp.push_back(singSuppR(i) - ifloor(singSuppR(i)));
    }
    temp.sort();
    temp.unique();
    
    flens::DenseVector<flens::Array<T> > singSupp(temp.size());
    FLENS_DEFAULT_INDEXTYPE i = 1;
    for (typename std::list<T>::const_iterator it = temp.begin(); it != temp.end(); ++it, ++i) {
        singSupp(i) = *it;
    }
    
    return singSupp;
}

template <typename T>
T
BSpline<T,Primal,Periodic,CDF>::tic(FLENS_DEFAULT_INDEXTYPE j) const
{
    return pow2i<T>(-j);
}

template <typename T>
flens::DenseVector<flens::Array<long double> > *
BSpline<T,Primal,Periodic,CDF>::getRefinement(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k, FLENS_DEFAULT_INDEXTYPE &refinement_j, FLENS_DEFAULT_INDEXTYPE &refinement_k_first,
											  FLENS_DEFAULT_INDEXTYPE &split, FLENS_DEFAULT_INDEXTYPE &refinement_k_restart) const
{
    refinement_j = j;

	refinement_k_restart = 1;

	//flens::DenseVector<flens::Array<long double> >* coeffs =  mra._RefCoeffs[0];
	// Left part
	if(k <= mra.rangeIL(j).lastIndex()){
		std::cerr << "Oops! I never expected there to be a left BSpline! " << std::endl;
		split = mra._periodicRefCoeffs[0].length() + 1;
		return &(mra._periodicRefCoeffs[0]);
	}
	// Inner part
	if(k <= mra.rangeII(j).lastIndex()){
        refinement_k_first = k + mra._periodicOffsets[0];
        split = mra._periodicRefCoeffs[0].length() + 1;
        return &(mra._periodicRefCoeffs[0]);
	}
	// Right part
	//flens::DenseVector<flens::Array<long double> >* coeffs =  mra._RefCoeffs[0];
	refinement_k_first = k + mra._rightOffsets[0];
	split = mra._split[k - mra.rangeIR(j).firstIndex()];
	refinement_k_restart = 1;
	return &(mra._rightRefCoeffs[0]);

}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
BSpline<T,Primal,Periodic,CDF>::getRefinementLevel(FLENS_DEFAULT_INDEXTYPE j) const
{
	return j;
}

template <typename T>
T
BSpline<T,Primal,Periodic,CDF>::getL2Norm(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    return mra._periodicL2Norms[0];
}

template <typename T>
T
BSpline<T,Primal,Periodic,CDF>::getH1SemiNorm(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    T pow2ij = (T)((FLENS_DEFAULT_INDEXTYPE) 1 << j);
    return pow2ij*mra._periodicH1SemiNorms[0];
}

template <typename T>
const flens::DenseVector<flens::Array<T> > &
BSpline<T,Primal,Periodic,CDF>::mask() const
{
    return phiR.a;
}

} // namespace lawa

