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

namespace lawa {

/*template <typename T>
Wavelet<T,Primal,Periodic,CDF>::Wavelet(int _d, int _d_)
    : d(_d), d_(_d_), mu(d&1),
      vanishingMoments(_d_), psiR(_d, _d_)
{
    assert(d<=d_);
    assert(((d+d_)&1)==0);
}

template <typename T>
Wavelet<T,Primal,Periodic,CDF>::Wavelet(const BSpline<T,Primal,Periodic,CDF> &_phi,
                                        const BSpline<T,Dual,Periodic,CDF> &_phi_)
    : d(_phi.d), d_(_phi_.d_), mu(d&1),
      vanishingMoments(d_), psiR(d, d_)
{
}
*/

template <typename T>
Wavelet<T,Primal,Periodic,CDF>::Wavelet(const Basis<T,Primal,Periodic,CDF> &_basis)
    : d(_basis.d), d_(_basis.d_), mu(d&1),
      vanishingMoments(d_), psiR(d,d_), basis(_basis)
{
}

template <typename T>
T
Wavelet<T,Primal,Periodic,CDF>::operator()(T x, int j, long k, unsigned short deriv) const
{
    // maximal support: [0,1]
    if((x < 0.) || (x > 1.)){
        return 0.;
    }
    
    // sum contributions of original spline on R
    // = 'wrapping' around [0,1]
    T val = 0;
    for(int l = ifloor(psiR.support(j,k).l1); l < iceil<int>(psiR.support(j,k).l2); ++l){
        val += psiR(l+x, j, k, deriv);
    }
    return val;
}

template <typename T>
PeriodicSupport<T>
Wavelet<T,Primal,Periodic,CDF>::support(int j, long k) const
{    
    Support<T> suppR = psiR.support(j,k);
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
Wavelet<T,Primal,Periodic,CDF>::singularSupport(int j, long k) const
{    
    if((psiR.support(j,k).l1 >= 0) && (psiR.support(j,k).l2 <= 1)){
         return linspace(support(j,k).l1, support(j,k).l2, 2*(d+d_)-1);
    }
    
    std::list<T> temp;
    flens::DenseVector<flens::Array<T> > singSuppR = linspace(psiR.support(j,k).l1, psiR.support(j,k).l2, 2*(d+d_)-1);
    temp.push_back(0.);
    temp.push_back(1.);
    for(int i = singSuppR.firstIndex(); i <= singSuppR.lastIndex(); ++i){
        temp.push_back(singSuppR(i) - ifloor(singSuppR(i)));
    }
    temp.sort();
    temp.unique();
    
    flens::DenseVector<flens::Array<T> > singSupp(temp.size());
    int i = 1;
    for (typename std::list<T>::const_iterator it = temp.begin(); it != temp.end(); ++it, ++i) {
        singSupp(i) = *it;
    }
    
    return singSupp;
}

template <typename T>
T
Wavelet<T,Primal,Periodic,CDF>::tic(int j) const
{
    return pow2i<T>(-(j+1));
}

template <typename T>
flens::DenseVector<flens::Array<long double> > *
Wavelet<T,Primal,Periodic,CDF>::getRefinement(int j, long k, int &refinement_j, long &refinement_k_first,
												long &split, long &refinement_k_restart) const
{
    refinement_j = j + 1;

	refinement_k_restart = 1;

	//flens::DenseVector<flens::Array<long double> >* coeffs =  mra._RefCoeffs[0];
	// Left part
	if(k <= basis.rangeJL(j).lastIndex()){
		std::cerr << "Ooops! Not implemented yet! " << std::endl;
    	split = basis._periodicRefCoeffs[0].length() + 1;
    	return &(basis._periodicRefCoeffs[0]);
	}
	// Inner part
	if(k <= basis.rangeJI(j).lastIndex()){
        refinement_k_first = 2*k + basis._innerOffsets[0];
    	split = basis._periodicRefCoeffs[0].length() + 1;
    	return &(basis._periodicRefCoeffs[0]);
	}
	// Right part
	//flens::DenseVector<flens::Array<long double> >* coeffs =  mra._RefCoeffs[0];
	int type = k - basis.rangeJR(j).firstIndex();
	refinement_k_first = 2*k + basis._innerOffsets[0];
	split = basis._split[type];
	refinement_k_restart = 1;
	return &(basis._rightRefCoeffs[type]);
}

template <typename T>
int
Wavelet<T,Primal,Periodic,CDF>::getRefinementLevel(int j) const {
	return j+1;
}

template <typename T>
T
Wavelet<T,Primal,Periodic,CDF>::getL2Norm(int j, long k) const
{
    return basis._periodicL2Norms[0];
}

template <typename T>
T
Wavelet<T,Primal,Periodic,CDF>::getH1SemiNorm(int j, long k) const
{
    T pow2ij = (T)(1L << j);
    return pow2ij*basis._periodicH1SemiNorms[0];
}

template <typename T>
const flens::DenseVector<flens::Array<T> > &
Wavelet<T,Primal,Periodic,CDF>::mask() const
{
    return psiR.b;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
Wavelet<T,Primal,Periodic,CDF>::mask(int d, int d_)
{   
    assert(d<=d_);
    assert(((d+d_)&1)==0);

    int mu = d & 1;
    BSpline<T,Dual,R,CDF> phi_(d,d_);
    flens::DenseVector<flens::Array<T> > b(flens::_(2-(d+mu)/2-d_, (d-mu)/2+d_));
    for (int k=b.firstIndex(); k<=b.lastIndex(); ++k) {
        int sign = (k&1) ? -1 : 1;
        b(k) = sign * phi_.a_(1-k);
    }
    return b;
}

} // namespace lawa

