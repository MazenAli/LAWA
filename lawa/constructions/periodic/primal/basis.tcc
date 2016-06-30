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

template <typename T>
Basis<T,Primal,Periodic,CDF>::Basis(FLENS_DEFAULT_INDEXTYPE _d, FLENS_DEFAULT_INDEXTYPE _d_, FLENS_DEFAULT_INDEXTYPE j)
    : d(_d), d_(_d_), j0(j), mra(d,d_,j), mra_(d,d_,j), 
      psi(*this), M1(psi), refinementbasis(_d, _d_, j), _j(j)
{
	if(d == 2 && d_ == 2){

        _periodicRefCoeffs = new flens::DenseVector<flens::Array<long double> >[1];
        _periodicRefCoeffs[0].engine().resize((FLENS_DEFAULT_INDEXTYPE)5,0);
        _periodicRefCoeffs[0] = 1.L/(4.L*std::sqrt(2.L)), 1.L/(2.L*std::sqrt(2.L)), - 3.L/(2.L*std::sqrt(2.L)), 1.L/(2.L*std::sqrt(2.L)), 1.L/(4.L*std::sqrt(2.L));

		_innerOffsets = new FLENS_DEFAULT_INDEXTYPE[1];
        _innerOffsets[0] = 0;

        _rightRefCoeffs = new flens::DenseVector<flens::Array<long double> >[2];
        _rightRefCoeffs[0].engine().resize((FLENS_DEFAULT_INDEXTYPE)6,0);
        _rightRefCoeffs[0] =  1.L/(4.L*std::sqrt(2.L)), 1.L/(2.L*std::sqrt(2.L)), - 3.L/(2.L*std::sqrt(2.L)), 1.L/(2.L*std::sqrt(2.L)), 1.L/(2.L*std::sqrt(2.L)), 1.L/(4.L*std::sqrt(2.L));

        _rightRefCoeffs[1].engine().resize((FLENS_DEFAULT_INDEXTYPE)6,0);
        _rightRefCoeffs[1] =  1.L/(4.L*std::sqrt(2.L)), 1.L/(2.L*std::sqrt(2.L)), 1.L/(2.L*std::sqrt(2.L)) , - 3.L/(2.L*std::sqrt(2.L)), 1.L/(2.L*std::sqrt(2.L)), 1.L/(4.L*std::sqrt(2.L));

		_split = new FLENS_DEFAULT_INDEXTYPE[2];
		_split[0] = 4;
		_split[1] = 2;

		_periodicL2Norms = new long double[1];
		_periodicL2Norms[0] = std::sqrt(3.L/4.L);
		_periodicH1SemiNorms = new long double[1];
		_periodicH1SemiNorms[0] = std::sqrt(16.5L);
	}
}

template <typename T>
Basis<T,Primal,Periodic, CDF>::~Basis()
{
	if(d == 2 && d_ == 2){
		delete[] _periodicRefCoeffs;
		delete[] _innerOffsets;
		delete[] _rightRefCoeffs;
		delete[] _split;
		delete[] _periodicL2Norms;
		delete[] _periodicH1SemiNorms;
	}
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
Basis<T,Primal,Periodic,CDF>::level() const
{
    return _j;
}

template <typename T>
void
Basis<T,Primal,Periodic,CDF>::setLevel(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=j0);
    _j = j;
}

template <typename T>
const BasisFunction<T,Primal,Periodic,CDF> &
Basis<T,Primal,Periodic,CDF>::generator(XType xtype) const
{
    if (xtype==XBSpline) {
        return mra.phi;
    } else {
        return psi;
    }
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
Basis<T,Primal,Periodic,CDF>::cardJ(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=j0);

    return pow2i<T>(j);
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
Basis<T,Primal,Periodic,CDF>::cardJL(FLENS_DEFAULT_INDEXTYPE /*j*/) const
{
	return std::max( std::ceil((d + d_ - 2)/2.0 - 1), 0.0);
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
Basis<T,Primal,Periodic,CDF>::cardJI(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=j0);

    return std::max(pow2i<T>(j) - cardJL() - cardJR(), 0.0);
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
Basis<T,Primal,Periodic,CDF>::cardJR(FLENS_DEFAULT_INDEXTYPE /*j*/) const
{
    return std::ceil((d + d_)/2.0 - 1) + 1;
}

template <typename T>
const flens::Range<FLENS_DEFAULT_INDEXTYPE>
Basis<T,Primal,Periodic,CDF>::rangeJ(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=j0);

    return flens::Range<FLENS_DEFAULT_INDEXTYPE>(1,pow2i<T>(j));
}


template <typename T>
const flens::Range<FLENS_DEFAULT_INDEXTYPE>
Basis<T,Primal,Periodic,CDF>::rangeJL(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=j0);

    return flens::Range<FLENS_DEFAULT_INDEXTYPE>(1,cardJL());
}


template <typename T>
const flens::Range<FLENS_DEFAULT_INDEXTYPE>
Basis<T,Primal,Periodic,CDF>::rangeJI(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=j0);

    return flens::Range<FLENS_DEFAULT_INDEXTYPE>(cardJL() + 1,pow2i<T>(j) - cardJR());
}

template <typename T>
const flens::Range<FLENS_DEFAULT_INDEXTYPE>
Basis<T,Primal,Periodic,CDF>::rangeJR(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=j0);

    return flens::Range<FLENS_DEFAULT_INDEXTYPE>(pow2i<T>(j) - std::ceil((d + d_)/2.0 - 1), pow2i<T>(j));
}

template <typename T>
void
Basis<T,Primal,Periodic,CDF>::getScalingNeighborsForScaling(FLENS_DEFAULT_INDEXTYPE j_scaling1, FLENS_DEFAULT_INDEXTYPE k_scaling1,
							  const Basis<T,Primal,Periodic,CDF> &/*secondbasis*/,
							  FLENS_DEFAULT_INDEXTYPE &j_scaling2, FLENS_DEFAULT_INDEXTYPE &k_scaling_first,
							  FLENS_DEFAULT_INDEXTYPE &k_scaling_last) const
{
    j_scaling2 = j_scaling1;

    FLENS_DEFAULT_INDEXTYPE kminR = k_scaling1 - d + 1;
    FLENS_DEFAULT_INDEXTYPE kmaxR = k_scaling1 + d - 1;

    k_scaling_first = kminR >= 1? kminR : (FLENS_DEFAULT_INDEXTYPE)mra.rangeI(j_scaling2).lastIndex() + kminR;
    k_scaling_last = kmaxR <= (FLENS_DEFAULT_INDEXTYPE)mra.rangeI(j_scaling2).lastIndex()? kmaxR : kmaxR - (FLENS_DEFAULT_INDEXTYPE)mra.rangeI(j_scaling2).lastIndex();

    if(k_scaling_first == k_scaling_last || k_scaling_first == k_scaling_last+1){
    	k_scaling_first = (FLENS_DEFAULT_INDEXTYPE)mra.rangeI(j_scaling2).firstIndex();
    	k_scaling_last = (FLENS_DEFAULT_INDEXTYPE)mra.rangeI(j_scaling2).lastIndex();
    }
}

template <typename T>
void
Basis<T,Primal,Periodic,CDF>::getScalingNeighborsForScaling(FLENS_DEFAULT_INDEXTYPE j_scaling1, FLENS_DEFAULT_INDEXTYPE k_scaling1,
							  const Basis<T,Primal,Interval,Dijkema> &secondbasis,
							  FLENS_DEFAULT_INDEXTYPE &j_scaling2, FLENS_DEFAULT_INDEXTYPE &k_scaling_first,
							  FLENS_DEFAULT_INDEXTYPE &k_scaling_last) const
{
    j_scaling2 = j_scaling1;

    // Use that k_interval = k_periodic + 1

    // Interval neighbours:
    FLENS_DEFAULT_INDEXTYPE k_int = k_scaling1 + 1;
    PeriodicSupport<T> supp = mra.phi.support(j_scaling1,k_scaling1);
    if(supp.gaplength() == 0){
        if (supp.l1==0.L) {
            k_scaling_first = (FLENS_DEFAULT_INDEXTYPE)secondbasis.mra.rangeI(j_scaling2).firstIndex();
            k_scaling_last  = std::min(k_int + d - 1, (FLENS_DEFAULT_INDEXTYPE)secondbasis.mra.rangeI(j_scaling2).lastIndex());
            return;
        }
        if (supp.l2<1.L) {
            k_scaling_first = std::max(k_int-d+1, (FLENS_DEFAULT_INDEXTYPE)secondbasis.mra.rangeI(j_scaling2).firstIndex());
            k_scaling_last  = std::min(k_int+d-1, (FLENS_DEFAULT_INDEXTYPE)secondbasis.mra.rangeI(j_scaling2).lastIndex());
            return;
        }
        k_scaling_first = std::max((FLENS_DEFAULT_INDEXTYPE)secondbasis.mra.rangeI(j_scaling2).firstIndex(),k_int - d + 1);
        k_scaling_last  = (FLENS_DEFAULT_INDEXTYPE)secondbasis.mra.rangeI(j_scaling2).lastIndex();
    }
    else{
    	// First Index: same as for 0 < l1 above
        k_scaling_first = std::max((FLENS_DEFAULT_INDEXTYPE)secondbasis.mra.rangeI(j_scaling2).firstIndex(),k_int - d + 1);
        // Last Index: similar as for l2 < 1 above, but for different k
        FLENS_DEFAULT_INDEXTYPE k_break = 1 + k_scaling1 - mra.rangeIR(j_scaling1).firstIndex();
        k_scaling_last = std::min(k_break+d-1, (FLENS_DEFAULT_INDEXTYPE)secondbasis.mra.rangeI(j_scaling2).lastIndex());

        if(k_scaling_first <= k_scaling_last+1){
        	k_scaling_first = (FLENS_DEFAULT_INDEXTYPE)secondbasis.mra.rangeI(j_scaling2).firstIndex();
        	k_scaling_last = (FLENS_DEFAULT_INDEXTYPE)secondbasis.mra.rangeI(j_scaling2).lastIndex();
        }
    }

    if(k_scaling_first == k_scaling_last || k_scaling_first == k_scaling_last+1){
    	k_scaling_first = (FLENS_DEFAULT_INDEXTYPE)secondbasis.mra.rangeI(j_scaling2).firstIndex();
    	k_scaling_last = (FLENS_DEFAULT_INDEXTYPE)secondbasis.mra.rangeI(j_scaling2).lastIndex();
    }
}

template <typename T>
void
Basis<T,Primal,Periodic,CDF>::getWaveletNeighborsForScaling(FLENS_DEFAULT_INDEXTYPE j_scaling1, FLENS_DEFAULT_INDEXTYPE k_scaling1,
							  const Basis<T,Primal,Periodic,CDF> &/*secondbasis*/,
							  FLENS_DEFAULT_INDEXTYPE &j_wavelet, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
							  FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const
{
    j_wavelet = j_scaling1;

    FLENS_DEFAULT_INDEXTYPE kminR = k_scaling1 - d + std::floor(0.5*((d&1)-d_)) + 1;
    FLENS_DEFAULT_INDEXTYPE kmaxR = k_scaling1 + d - 1 + std::ceil(0.5*((d&1)+d_)) - 1;

    k_wavelet_first = kminR >= 1? kminR : rangeJ(j_wavelet).lastIndex() + kminR;
    k_wavelet_last = kmaxR <= rangeJ(j_wavelet).lastIndex()? kmaxR : kmaxR - rangeJ(j_wavelet).lastIndex();

    if(k_wavelet_first == k_wavelet_last || k_wavelet_first == k_wavelet_last+1){
    	k_wavelet_first = rangeJ(j_wavelet).firstIndex();
    	k_wavelet_last = rangeJ(j_wavelet).lastIndex();
    }
}

template <typename T>
void
Basis<T,Primal,Periodic,CDF>::getWaveletNeighborsForScaling(FLENS_DEFAULT_INDEXTYPE j_scaling1, FLENS_DEFAULT_INDEXTYPE k_scaling1,
							  const Basis<T,Primal,Interval,Dijkema> &secondbasis,
							  FLENS_DEFAULT_INDEXTYPE &j_wavelet, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
							  FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const
{
    j_wavelet = j_scaling1;

    // Use that k_interval = k_periodic + 1
    // Interval neighbours:
    FLENS_DEFAULT_INDEXTYPE k_int = k_scaling1 + 1;

    PeriodicSupport<T> supp = mra.phi.support(j_scaling1,k_scaling1);

    if(supp.gaplength()==0){
        if (supp.l1==0.L) {
            k_wavelet_first = 1;
            k_wavelet_last = k_wavelet_first + secondbasis.cardJL(j_wavelet) + d/2;
            k_wavelet_last  = std::min(k_wavelet_last, (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJR(j_wavelet).lastIndex());
            return;
        }
        if (0<supp.l1 && supp.l2<1.L) {
            k_wavelet_first  = std::max((FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJL(j_wavelet).firstIndex(), k_int - (d+d_) + 1);
            k_wavelet_last   = std::min((FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJR(j_wavelet).lastIndex(),  k_int + (d+d_) - 1);
            return;
        }
        k_wavelet_last   = secondbasis.rangeJ(j_wavelet).lastIndex();
        k_wavelet_first  = k_wavelet_last - (secondbasis.cardJR(j_wavelet) + d/2) + 1;
        k_wavelet_first  = std::max((FLENS_DEFAULT_INDEXTYPE) 1, k_wavelet_first);
    }
    else{
        k_wavelet_first  = std::max((FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJL(j_wavelet).firstIndex(), k_int - (d+d_) + 1);
        FLENS_DEFAULT_INDEXTYPE k_break = 1 + k_scaling1 - mra.rangeIR(j_scaling1).firstIndex();
        k_wavelet_last   = std::min((FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJR(j_wavelet).lastIndex(),  k_break + (d+d_) - 1);

        if(k_wavelet_first <= k_wavelet_last+1){
        	k_wavelet_first = (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet).firstIndex();
        	k_wavelet_last = (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet).lastIndex();
        	return;
        }

    }

    if(k_wavelet_first == k_wavelet_last || k_wavelet_first == k_wavelet_last+1){
    	k_wavelet_first = secondbasis.rangeJ(j_wavelet).firstIndex();
    	k_wavelet_last = secondbasis.rangeJ(j_wavelet).lastIndex();
    }
}

template <typename T>
template <typename SecondRefinementBasis>
void
Basis<T,Primal,Periodic,CDF>::getBSplineNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet, FLENS_DEFAULT_INDEXTYPE k_wavelet,
                              const SecondRefinementBasis &secondrefinementbasis,
                              FLENS_DEFAULT_INDEXTYPE &j_bspline, FLENS_DEFAULT_INDEXTYPE &k_bspline_first,
                              FLENS_DEFAULT_INDEXTYPE &k_bspline_last) const
{
    ct_assert(SecondRefinementBasis::Side==Primal and SecondRefinementBasis::Domain==Interval and SecondRefinementBasis::Cons==Dijkema);

    j_bspline = j_wavelet;

    FLENS_DEFAULT_INDEXTYPE kminR = k_wavelet - d + 1 - std::ceil(0.5*((d&1)+d_)) + 1;
    FLENS_DEFAULT_INDEXTYPE kmaxR = k_wavelet + d + std::ceil(0.5*(d_-(d&1))) - 1;

    k_bspline_first = kminR >= 1? kminR : (FLENS_DEFAULT_INDEXTYPE)mra.rangeI(j_bspline).lastIndex() + kminR;
    k_bspline_last = kmaxR <= (FLENS_DEFAULT_INDEXTYPE)mra.rangeI(j_bspline).lastIndex()? kmaxR : kmaxR - (FLENS_DEFAULT_INDEXTYPE)mra.rangeI(j_bspline).lastIndex();

    /*
     * TODO: Extend to d > 2
     *
     * (Assumption: d == 2 !!!!)
     */


    // Take all Bsplines on this level? Then almost no further adjustements are necessary...
    if(k_bspline_first == k_bspline_last || k_bspline_first == k_bspline_last+1){

    	k_bspline_first = (FLENS_DEFAULT_INDEXTYPE)secondrefinementbasis.mra.rangeI(j_bspline).firstIndex();
    	k_bspline_last = (FLENS_DEFAULT_INDEXTYPE)secondrefinementbasis.mra.rangeI(j_bspline).lastIndex();

    }
    else{
        // If not, we have to shift the translation indizes
    	FLENS_DEFAULT_INDEXTYPE firstIndex = (FLENS_DEFAULT_INDEXTYPE)secondrefinementbasis.mra.rangeI(j_bspline).firstIndex();
    	FLENS_DEFAULT_INDEXTYPE lastIndex = (FLENS_DEFAULT_INDEXTYPE)secondrefinementbasis.mra.rangeI(j_bspline).lastIndex();
        k_bspline_first = std::max(std::min(k_bspline_first+1, lastIndex), firstIndex);
        k_bspline_last  = std::max(std::min(k_bspline_last+1, lastIndex), firstIndex);
    }

    // Test if we can drop the left or right boundary spline
    PeriodicSupport<T> supp_w = psi.support(j_wavelet, k_wavelet);

    if(overlap(supp_w, secondrefinementbasis.mra.phi.support(j_bspline, k_bspline_first)) <= 0){
    	k_bspline_first = k_bspline_first == secondrefinementbasis.mra.rangeI(j_bspline).lastIndex()? secondrefinementbasis.mra.rangeI(j_bspline).firstIndex() : k_bspline_first+1;
    }
    if(overlap(supp_w, secondrefinementbasis.mra.phi.support(j_bspline, k_bspline_last)) <= 0){
    	k_bspline_last = k_bspline_last == secondrefinementbasis.mra.rangeI(j_bspline).firstIndex()? secondrefinementbasis.mra.rangeI(j_bspline).lastIndex() : k_bspline_last-1;
    }
}


template <typename T>
template <typename SecondBasis>
void
Basis<T,Primal,Periodic,CDF>::getScalingNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet, FLENS_DEFAULT_INDEXTYPE k_wavelet,
                              const SecondBasis &secondbasis,
                              FLENS_DEFAULT_INDEXTYPE &j_scaling, FLENS_DEFAULT_INDEXTYPE &k_scaling_first,
                              FLENS_DEFAULT_INDEXTYPE &k_scaling_last) const
{
    ct_assert(SecondBasis::Side==Primal and SecondBasis::Domain==Periodic and SecondBasis::Cons==CDF);
    j_scaling = j_wavelet;

    FLENS_DEFAULT_INDEXTYPE kminR = k_wavelet - d + 1 - std::ceil(0.5*((d&1)+d_)) + 1;
    FLENS_DEFAULT_INDEXTYPE kmaxR = k_wavelet + d + std::ceil(0.5*(d_-(d&1))) - 1;

    k_scaling_first = kminR >= 1? kminR : (FLENS_DEFAULT_INDEXTYPE)mra.rangeI(j_scaling).lastIndex() + kminR;
    k_scaling_last = kmaxR <= (FLENS_DEFAULT_INDEXTYPE)mra.rangeI(j_scaling).lastIndex()? kmaxR : kmaxR - (FLENS_DEFAULT_INDEXTYPE)mra.rangeI(j_scaling).lastIndex();

    if(k_scaling_first == k_scaling_last || k_scaling_first == k_scaling_last+1){
    	k_scaling_first = (FLENS_DEFAULT_INDEXTYPE)secondbasis.mra.rangeI(j_scaling).firstIndex();
    	k_scaling_last = (FLENS_DEFAULT_INDEXTYPE)secondbasis.mra.rangeI(j_scaling).lastIndex();
    }
}


template <typename T>
void
Basis<T,Primal,Periodic,CDF>::getWaveletNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet1, FLENS_DEFAULT_INDEXTYPE k_wavelet1,
                              const Basis<T,Primal,Periodic,CDF> &/*secondbasis*/,
                              FLENS_DEFAULT_INDEXTYPE &j_wavelet2, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                              FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const
{
    j_wavelet2 = j_wavelet1;


    FLENS_DEFAULT_INDEXTYPE kminR = k_wavelet1 - d - d_ + 1 + 1;
    FLENS_DEFAULT_INDEXTYPE kmaxR = k_wavelet1 + d + d_ - 1 - 1;

    k_wavelet_first = kminR >= 1? kminR : rangeJ(j_wavelet2).lastIndex() + kminR;
    k_wavelet_last = kmaxR <= rangeJ(j_wavelet2).lastIndex()? kmaxR : kmaxR - rangeJ(j_wavelet2).lastIndex();

    if(k_wavelet_first == k_wavelet_last || k_wavelet_first == k_wavelet_last+1){
    	k_wavelet_first = rangeJ(j_wavelet2).firstIndex();
    	k_wavelet_last = rangeJ(j_wavelet2).lastIndex();
    }
}

template <typename T>
void
Basis<T,Primal,Periodic,CDF>::getWaveletNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet1, FLENS_DEFAULT_INDEXTYPE k_wavelet1,
                              const Basis<T,Primal,Interval,Dijkema> &secondbasis,
                              FLENS_DEFAULT_INDEXTYPE &j_wavelet2, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                              FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const
{
    j_wavelet2 = j_wavelet1;

    PeriodicSupport<T> supp = psi.support(j_wavelet1,k_wavelet1);
    if(supp.gaplength()==0){
        if (supp.l1==0.L) {
            k_wavelet_first = (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).firstIndex();
            k_wavelet_last  = std::min(k_wavelet1+2*d, (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).lastIndex());
            return;
        }
        if (supp.l2<1.L) {
            k_wavelet_first = std::max(k_wavelet1-2*d+1, (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).firstIndex());
            k_wavelet_last  = std::min(k_wavelet1+2*d-1, (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).lastIndex());
            return;
        }
        k_wavelet_first = std::max((FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).firstIndex(),k_wavelet1 - 2*d);
        k_wavelet_last  = (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).lastIndex();
    }
    else{
        k_wavelet_first = std::max((FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).firstIndex(),k_wavelet1+1-2*d+1);
        k_wavelet_last  = std::min((FLENS_DEFAULT_INDEXTYPE)2*d-1, (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).lastIndex());

        if(k_wavelet_first <= k_wavelet_last+1){
        	k_wavelet_first = (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).firstIndex();
        	k_wavelet_last = (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).lastIndex();
        	return;
        }
    }

    if(k_wavelet_first == k_wavelet_last || k_wavelet_first == k_wavelet_last+1){
    	k_wavelet_first = secondbasis.rangeJ(j_wavelet2).firstIndex();
    	k_wavelet_last = secondbasis.rangeJ(j_wavelet2).lastIndex();
    }
}

template <typename T>
template <typename SecondBasis>
void
Basis<T,Primal,Periodic,CDF>::getLowerWaveletNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet1, FLENS_DEFAULT_INDEXTYPE k_wavelet1,
                                       const SecondBasis &/*secondbasis*/,
                                       FLENS_DEFAULT_INDEXTYPE &j_wavelet2, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                       FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const
{
    ct_assert(SecondBasis::Side==Primal and SecondBasis::Domain==Periodic and SecondBasis::Cons==CDF);
    j_wavelet2 = j_wavelet1-1;


    FLENS_DEFAULT_INDEXTYPE kminR = std::floor(0.5*k_wavelet1 + 0.5 - 0.75*(d + d_)) + 1;
    FLENS_DEFAULT_INDEXTYPE kmaxR = std::ceil(0.5*k_wavelet1 - 1 + 0.75*(d+d_)) - 1;

    k_wavelet_first = kminR >= 1? kminR : rangeJ(j_wavelet2).lastIndex() + kminR;
    k_wavelet_last = kmaxR <= rangeJ(j_wavelet2).lastIndex()? kmaxR : kmaxR - rangeJ(j_wavelet2).lastIndex();

    if(k_wavelet_first == k_wavelet_last || k_wavelet_first == k_wavelet_last+1){
    	k_wavelet_first = rangeJ(j_wavelet2).firstIndex();
    	k_wavelet_last = rangeJ(j_wavelet2).lastIndex();
    }
}

template <typename T>
void
Basis<T,Primal,Periodic,CDF>::getHigherWaveletNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet1, FLENS_DEFAULT_INDEXTYPE k_wavelet1,
                                   const Basis<T,Primal,Periodic,CDF> &/*secondbasis*/,
                                   FLENS_DEFAULT_INDEXTYPE &j_wavelet2, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                   FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const
{
    j_wavelet2 = j_wavelet1+1;


    FLENS_DEFAULT_INDEXTYPE kminR = 2 + 2*k_wavelet1 - std::ceil(1.5*(d + d_)) + 1;
    FLENS_DEFAULT_INDEXTYPE kmaxR = 2*k_wavelet1 - 1 + std::ceil(1.5*(d+d_)) - 1;

    k_wavelet_first = kminR >= 1? kminR : rangeJ(j_wavelet2).lastIndex() + kminR;
    k_wavelet_last = kmaxR <= rangeJ(j_wavelet2).lastIndex()? kmaxR : kmaxR - rangeJ(j_wavelet2).lastIndex();

    if(k_wavelet_first == k_wavelet_last || k_wavelet_first == k_wavelet_last+1){
    	k_wavelet_first = rangeJ(j_wavelet2).firstIndex();
    	k_wavelet_last = rangeJ(j_wavelet2).lastIndex();
    }
}

template <typename T>
void
Basis<T,Primal,Periodic,CDF>::getHigherWaveletNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet1, FLENS_DEFAULT_INDEXTYPE k_wavelet1,
                                   const Basis<T,Primal,Interval,Dijkema> &secondbasis,
                                   FLENS_DEFAULT_INDEXTYPE &j_wavelet2, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                   FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const
{
    j_wavelet2 = j_wavelet1+1;

    // Use that k_interval = k_periodic + 1
    FLENS_DEFAULT_INDEXTYPE k_tilde = (k_wavelet1 + 1)*2;

    // Interval calculations
    PeriodicSupport<T> supp = psi.support(j_wavelet1,k_wavelet1);
    if(supp.gaplength() == 0){
        if (supp.l1==0.L) {
            k_wavelet_first = (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).firstIndex();
            k_wavelet_last  = std::min(k_tilde + 3*d, (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).lastIndex());
            return;
        }
        if (supp.l2<1.L) {
            k_wavelet_first = std::max(k_tilde-2*d-1, (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).firstIndex());
            k_wavelet_last  = std::min(k_tilde+2*d+1, (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).lastIndex());
            return;
        }
        k_wavelet_first = std::max((FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).firstIndex(),k_tilde - 3*d);
        k_wavelet_last  = (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).lastIndex();
    }
    else{
        k_wavelet_first = std::max((FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).firstIndex(),k_tilde - 3*d);
        FLENS_DEFAULT_INDEXTYPE k_break = 1 + k_wavelet1 - rangeJR(j_wavelet1).firstIndex();
        k_wavelet_last  = std::min(k_break+3*d, (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).lastIndex());

        if(k_wavelet_first <= k_wavelet_last+1){
        	k_wavelet_first = (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).firstIndex();
        	k_wavelet_last = (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).lastIndex();
        	return;
        }
    }

    if(k_wavelet_first == k_wavelet_last || k_wavelet_first == k_wavelet_last+1){
    	k_wavelet_first = (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).firstIndex();
    	k_wavelet_last = (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).lastIndex();
    }
}

} // namespace lawa

