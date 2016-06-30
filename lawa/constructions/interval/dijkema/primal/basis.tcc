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
#include <lawa/aux/compiletime_assert.h>
#include <lawa/constructions/interval/initial_stable_completion.h>

namespace lawa {

template <typename T>
Basis<T,Primal,Interval,Dijkema>::Basis(FLENS_DEFAULT_INDEXTYPE _d, FLENS_DEFAULT_INDEXTYPE _d_, FLENS_DEFAULT_INDEXTYPE j)
    : mra(_d, j), mra_(_d, _d_, j),
      d(_d), d_(_d_), mu(d&1),
      min_j0(mra_.min_j0), j0(mra_.j0), _bc(2,0), _j(-1), psi(*this), refinementbasis(*this),
      LaplaceOp1D(_d,*this), IdentityOp1D(_d,*this)
{
    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Mj1, Mj1_;
    initial_stable_completion(mra,mra_,Mj1,Mj1_);
    const FLENS_DEFAULT_INDEXTYPE cons_j = ((d==2) && ((d_==2)||(d_==4))) ? mra_.min_j0+1 : mra_.min_j0;
    M1 = flens::RefinementMatrix<T,Interval,Dijkema>(d+d_-2, d+d_-2, Mj1, min_j0, cons_j);
    _j = std::max(min_j0,j);
    setLevel(_j);

    if (d==2 && d_==2) {
        // left part
        _leftRefCoeffs = new flens::DenseVector<flens::Array<long double> >[2];
        _leftRefCoeffs[0].engine().resize((FLENS_DEFAULT_INDEXTYPE)4,0);
        _leftRefCoeffs[0] =  -1.060660171779819194, 0.795495128834866838, -0.176776695296636810, -0.088388347648318405;
        _leftRefCoeffs[1].engine().resize((FLENS_DEFAULT_INDEXTYPE)5,0);
        _leftRefCoeffs[1] =  -0.1767766952966368, -0.3535533905932737, 1.0606601717798211, -0.3535533905932737, -0.1767766952966368;
        _leftOffsets = new FLENS_DEFAULT_INDEXTYPE[2];
        _leftOffsets[0] =  1;
        _leftOffsets[1] =  2;
        _leftL2Norms = new long double[2];
        _leftL2Norms[0] =  std::sqrt(0.5L);   _leftL2Norms[1] =  std::sqrt(0.75L);
        _leftH1SemiNorms  = new long double[2];
        _leftH1SemiNorms[0] =  std::sqrt(17.625L);
        _leftH1SemiNorms[1] =  std::sqrt(16.5L);

        // inner part
        _innerRefCoeffs = new flens::DenseVector<flens::Array<long double> >[1];
        _innerRefCoeffs[0].engine().resize((FLENS_DEFAULT_INDEXTYPE)5,0);
        _innerRefCoeffs[0] =  -0.1767766952966368, -0.3535533905932737, 1.0606601717798211, -0.3535533905932737, -0.1767766952966368;
        _innerOffsets = new FLENS_DEFAULT_INDEXTYPE[1];
        _innerOffsets[0] =  -2;
        _innerL2Norms = new long double[1];
        _innerL2Norms[0] =  std::sqrt(0.75L);
        _innerH1SemiNorms  = new long double[1];
        _innerH1SemiNorms[0] =  std::sqrt(16.5L);

        // inner part
        _rightRefCoeffs = new flens::DenseVector<flens::Array<long double> >[2];
        _rightRefCoeffs[0].engine().resize((FLENS_DEFAULT_INDEXTYPE)5,0);
        _rightRefCoeffs[0] =  -0.1767766952966368, -0.3535533905932737, 1.0606601717798211, -0.3535533905932737, -0.1767766952966368;
        _rightRefCoeffs[1].engine().resize((FLENS_DEFAULT_INDEXTYPE)4,0);
        _rightRefCoeffs[1] =  -0.088388347648318613, -0.176776695296637226, 0.795495128834866505, -1.060660171779819638;
        _rightOffsets = new FLENS_DEFAULT_INDEXTYPE[2];
        _rightOffsets[0] =  - 2;
        _rightOffsets[1] =    0;
        _rightL2Norms = new long double[2];
        _rightL2Norms[0] =  std::sqrt(0.75L);   _rightL2Norms[1] =  std::sqrt(0.5L);
        _rightH1SemiNorms  = new long double[2];
        _rightH1SemiNorms[0] =  std::sqrt(16.5L);
        _rightH1SemiNorms[1] =  std::sqrt(17.625L);
    }
}

template <typename T>
Basis<T,Primal,Interval,Dijkema>::~Basis()
{
    if ((d==2 && d_==2) || (_bc(0)==1 && _bc(1)==1 && d==3 && d_==3)) {
        delete[] _leftRefCoeffs,
        delete[] _innerRefCoeffs,
        delete[] _rightRefCoeffs;
        delete[] _leftOffsets,
        delete[] _innerOffsets,
        delete[] _rightOffsets;
        delete[] _leftL2Norms;
        delete[] _innerL2Norms;
        delete[] _rightL2Norms;
        delete[] _leftH1SemiNorms;
        delete[] _innerH1SemiNorms;
        delete[] _rightH1SemiNorms;
    }

}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
Basis<T,Primal,Interval,Dijkema>::level() const
{
    return _j;
}

template <typename T>
void
Basis<T,Primal,Interval,Dijkema>::setLevel(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=min_j0);
    _j = j;
    M1.setLevel(_j);
    mra.setLevel(_j);
    mra_.setLevel(_j);
}

template <typename T>
template <BoundaryCondition BC>
void
Basis<T,Primal,Interval,Dijkema>::enforceBoundaryCondition()
{
    if ((_bc(0)==0) && (_bc(1)==0)) {
        _bc(0) = _bc(1) = 1;
        mra.template enforceBoundaryCondition<BC>();
        mra_.template enforceBoundaryCondition<BC>();
        flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > Mj1, Mj1_;
        initial_stable_completion(mra,mra_,Mj1,Mj1_);
        const FLENS_DEFAULT_INDEXTYPE cons_j = ((d==2) && ((d_==2)||(d_==4))) ? mra_.min_j0+1 : mra_.min_j0;
        M1.left.resize(0);
        M1.right.resize(0);
        M1.leftband.resize(0);
        M1.rightband.resize(0);
        M1.lengths.resize(0);

        M1 = flens::RefinementMatrix<T,Interval,Dijkema>(d+d_-2, d+d_-2, Mj1, min_j0, cons_j);
        setLevel(_j);
    }
    // Refinement coefficients only in double prec. available due to missing support of higher
    // precision in blas routines.
    if (d==2 && d_==2) {

        delete[] _leftRefCoeffs,
        delete[] _innerRefCoeffs,
        delete[] _rightRefCoeffs;
        delete[] _leftOffsets,
        delete[] _innerOffsets,
        delete[] _rightOffsets;
        delete[] _leftL2Norms;
        delete[] _innerL2Norms;
        delete[] _rightL2Norms;
        delete[] _leftH1SemiNorms;
        delete[] _innerH1SemiNorms;
        delete[] _rightH1SemiNorms;

        // left part
        _leftRefCoeffs = new flens::DenseVector<flens::Array<long double> >[2];
        _leftRefCoeffs[0].engine().resize((FLENS_DEFAULT_INDEXTYPE)3,0);
        _leftRefCoeffs[0] =  1.2374368670764579, -0.3535533905932737, -0.1767766952966368;
        _leftRefCoeffs[1].engine().resize((FLENS_DEFAULT_INDEXTYPE)5,0);
        _leftRefCoeffs[1] =  -0.1767766952966368, -0.3535533905932737, 1.0606601717798211, -0.3535533905932737, -0.1767766952966368;
        _leftOffsets = new FLENS_DEFAULT_INDEXTYPE[2];
        _leftOffsets[0] =  2;
        _leftOffsets[1] =  2;
        _leftL2Norms = new long double[2];
        _leftL2Norms[0] =  1.L;   _leftL2Norms[1] =  std::sqrt(0.75L);
        _leftH1SemiNorms  = new long double[2];
        _leftH1SemiNorms[0] =  std::sqrt(16.5L);   _leftH1SemiNorms[1] =  std::sqrt(16.5L);

        // inner part
        _innerRefCoeffs = new flens::DenseVector<flens::Array<long double> >[1];
        _innerRefCoeffs[0].engine().resize((FLENS_DEFAULT_INDEXTYPE)5,0);
        _innerRefCoeffs[0] =  -0.1767766952966368, -0.3535533905932737, 1.0606601717798211, -0.3535533905932737, -0.1767766952966368;
        _innerOffsets = new FLENS_DEFAULT_INDEXTYPE[1];
        _innerOffsets[0] =  -2;
        _innerL2Norms = new long double[1];
        _innerL2Norms[0] =  std::sqrt(0.75L);
        _innerH1SemiNorms  = new long double[1];
        _innerH1SemiNorms[0] =  std::sqrt(16.5L);

        // inner part
        _rightRefCoeffs = new flens::DenseVector<flens::Array<long double> >[2];
        _rightRefCoeffs[0].engine().resize((FLENS_DEFAULT_INDEXTYPE)5,0);
        _rightRefCoeffs[0] =  -0.1767766952966368, -0.3535533905932737, 1.0606601717798211, -0.3535533905932737, -0.1767766952966368;
        _rightRefCoeffs[1].engine().resize((FLENS_DEFAULT_INDEXTYPE)3,0);
        _rightRefCoeffs[1] =  -0.1767766952966368, -0.3535533905932737, 1.2374368670764579;
        _rightOffsets = new FLENS_DEFAULT_INDEXTYPE[2];
        _rightOffsets[0] =  - 2;
        _rightOffsets[1] =    0;
        _rightL2Norms = new long double[2];
        _rightL2Norms[0] =  std::sqrt(0.75L);   _rightL2Norms[1] =  1.L;
        _rightH1SemiNorms  = new long double[2];
        _rightH1SemiNorms[0] =  std::sqrt(16.5L); _rightH1SemiNorms[1] =  std::sqrt(16.5L);
    }
    else if (d==3 && d_==3) {

        // left part
        _leftRefCoeffs = new flens::DenseVector<flens::Array<long double> >[4];
        _leftRefCoeffs[0].engine().resize((FLENS_DEFAULT_INDEXTYPE)7,0);
        _leftRefCoeffs[0] =  0.825002583178343, -0.31011594165461, -0.063691491010539, 0.044358015602397, 0.014032578184197, -0.000847605393676, -0.0002825351312254;
        _leftRefCoeffs[1].engine().resize((FLENS_DEFAULT_INDEXTYPE)7,0);
        _leftRefCoeffs[1] = -0.084857408555486,  0.12118335399876,  0.702979696218226757, -0.707687538747675, -0.110818350898946, 0.1407121822690638, 0.0469040607563546;
        _leftRefCoeffs[2].engine().resize((FLENS_DEFAULT_INDEXTYPE)8,0);
        _leftRefCoeffs[2] = -0.046875, -0.140625, 0.109375, 0.703125, -0.703125, -0.109375, 0.140625, 0.046875;
        _leftRefCoeffs[3].engine().resize((FLENS_DEFAULT_INDEXTYPE)8,0);
        _leftRefCoeffs[3] = -0.046875, -0.140625, 0.109375, 0.703125, -0.703125, -0.109375, 0.140625, 0.046875;
        _leftOffsets = new FLENS_DEFAULT_INDEXTYPE[4];
        _leftOffsets[0] =  2;
        _leftOffsets[1] =  2;
        _leftOffsets[2] =  3;
        _leftOffsets[3] =  5;
        _leftL2Norms = new long double[4];
        _leftL2Norms[0] =  std::sqrt(0.18306923631658L);
        _leftL2Norms[1] =  std::sqrt(0.4182530928956412L);
        _leftL2Norms[2] =  std::sqrt(0.419921875L);
        _leftL2Norms[3] =  std::sqrt(0.419921875L);
        _leftH1SemiNorms  = new long double[4];
        _leftH1SemiNorms[0] =  std::sqrt(283.1110396147645929L/64.L);
        _leftH1SemiNorms[1] =  std::sqrt(357.10493898574242166L/64.L);
        _leftH1SemiNorms[2] =  std::sqrt(362.5L/64.L);
        _leftH1SemiNorms[3] =  std::sqrt(362.5L/64.L);


        // inner part
        _innerRefCoeffs = new flens::DenseVector<flens::Array<long double> >[2];
        _innerRefCoeffs[0].engine().resize((FLENS_DEFAULT_INDEXTYPE)8,0);
        _innerRefCoeffs[0] =  0.046875,  0.140625, -0.109375, -0.703125, 0.703125, 0.109375, -0.140625, -0.046875;
        _innerRefCoeffs[1].engine().resize((FLENS_DEFAULT_INDEXTYPE)8,0);
        _innerRefCoeffs[1] = -0.046875, -0.140625, 0.109375, 0.703125, -0.703125, -0.109375, 0.140625, 0.046875;
        _innerOffsets = new FLENS_DEFAULT_INDEXTYPE[2];
        _innerOffsets[0] =  -3;
        _innerOffsets[1] =  -3;
        _innerL2Norms = new long double[1];
        _innerL2Norms[0] =  std::sqrt(0.419921875L);
        _innerH1SemiNorms  = new long double[1];
        _innerH1SemiNorms[0] =  std::sqrt(362.5L/64.L);

        // inner part
        _rightRefCoeffs = new flens::DenseVector<flens::Array<long double> >[4];
        _rightRefCoeffs[0].engine().resize((FLENS_DEFAULT_INDEXTYPE)8,0);
        _rightRefCoeffs[0] =  0.046875, 0.140625, -0.109375, -0.703125, 0.703125, 0.109375, -0.140625, -0.046875;
        _rightRefCoeffs[1].engine().resize((FLENS_DEFAULT_INDEXTYPE)8,0);
        _rightRefCoeffs[1] =  0.046875, 0.140625, -0.109375, -0.703125, 0.703125, 0.109375, -0.140625, -0.046875;
        _rightRefCoeffs[2].engine().resize((FLENS_DEFAULT_INDEXTYPE)7,0);
        _rightRefCoeffs[2] = 0.0469040607563546, 0.1407121822690638,  -0.110818350898946,  -0.707687538747675, 0.702979696218226757, 0.12118335399876, -0.084857408555486;
        _rightRefCoeffs[3].engine().resize((FLENS_DEFAULT_INDEXTYPE)7,0);
        _rightRefCoeffs[3] = -0.0002825351312254, -0.000847605393676, 0.014032578184197, 0.044358015602397, -0.063691491010539, -0.31011594165461, 0.825002583178343;
        _rightOffsets = new FLENS_DEFAULT_INDEXTYPE[4];
        _rightOffsets[0] =  - 7;
        _rightOffsets[1] =  - 5;
        _rightOffsets[2] =  - 3;
        _rightOffsets[3] =  - 3;
        _rightL2Norms = new long double[4];
        _rightL2Norms[0] =  std::sqrt(0.419921875L);
        _rightL2Norms[1] =  std::sqrt(0.419921875L);
        _rightL2Norms[2] =  std::sqrt(0.4182530928956412L);
        _rightL2Norms[3] =  std::sqrt(0.18306923631658L);
        _rightH1SemiNorms  = new long double[4];
        _rightH1SemiNorms[0] =  std::sqrt(362.5L/64.L);
        _rightH1SemiNorms[1] =  std::sqrt(362.5L/64.L);
        _rightH1SemiNorms[2] =  std::sqrt(357.10493898574242166L/64.L);
        _rightH1SemiNorms[3] =  std::sqrt(283.1110396147645929L/64.L);
    }
}

template <typename T>
const BasisFunction<T,Primal,Interval,Dijkema> &
Basis<T,Primal,Interval,Dijkema>::generator(XType xtype) const
{
    if (xtype==XBSpline) {
        return mra.phi; 
    } else {
        return psi;
    }
}

// cardinalities of whole, left, inner, right index sets (primal).
template <typename T>
FLENS_DEFAULT_INDEXTYPE
Basis<T,Primal,Interval,Dijkema>::cardJ(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=min_j0);
    return pow2i<T>(j);
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
Basis<T,Primal,Interval,Dijkema>::cardJL(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=min_j0);
    return d + d_ - 2;
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
Basis<T,Primal,Interval,Dijkema>::cardJI(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=min_j0);
    return pow2i<T>(j) - 2*(d + d_ - 2);
}

template <typename T>
FLENS_DEFAULT_INDEXTYPE
Basis<T,Primal,Interval,Dijkema>::cardJR(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=min_j0);
    return d + d_ - 2;
}

// ranges of whole, left, inner, right index sets (primal).
template <typename T>
const flens::Range<FLENS_DEFAULT_INDEXTYPE>
Basis<T,Primal,Interval,Dijkema>::rangeJ(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=min_j0);
    return flens::_(1,pow2i<T>(j));
}

template <typename T>
const flens::Range<FLENS_DEFAULT_INDEXTYPE>
Basis<T,Primal,Interval,Dijkema>::rangeJL(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=min_j0);
    return flens::_(1,d+d_-2);
}

template <typename T>
const flens::Range<FLENS_DEFAULT_INDEXTYPE>
Basis<T,Primal,Interval,Dijkema>::rangeJI(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=min_j0);
    return flens::_(d+d_-1,pow2i<T>(j)-(d+d_-2));
}

template <typename T>
const flens::Range<FLENS_DEFAULT_INDEXTYPE>
Basis<T,Primal,Interval,Dijkema>::rangeJR(FLENS_DEFAULT_INDEXTYPE j) const
{
    assert(j>=min_j0);
    return flens::_(pow2i<T>(j)-(d+d_-3),pow2i<T>(j));
}


template <typename T>
template <typename SecondBasis>
void
Basis<T,Primal,Interval,Dijkema>::
getWaveletNeighborsForBSpline(FLENS_DEFAULT_INDEXTYPE j_bspline, FLENS_DEFAULT_INDEXTYPE k_bspline,
                              const SecondBasis &secondbasis,
                              FLENS_DEFAULT_INDEXTYPE &j_wavelet, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first, FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const
{
    ct_assert(SecondBasis::Side==Primal and SecondBasis::Domain==Interval
              and SecondBasis::Cons==Dijkema);
    j_wavelet = j_bspline;
    Support<T> supp = mra.phi.support(j_bspline,k_bspline);
    T a = supp.l1, b = supp.l2;

    if (a==0.L) {
        k_wavelet_first = 1;
        k_wavelet_last = k_wavelet_first + secondbasis.cardJL(j_wavelet) + d/2;
        k_wavelet_last  = std::min(k_wavelet_last, (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJR(j_wavelet).lastIndex());
        return;
    }
    if (0<a && b<1.L) {
        k_wavelet_first  = std::max((FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJL(j_wavelet).firstIndex(), k_bspline - (d+d_) + 1);
        k_wavelet_last   = std::min((FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJR(j_wavelet).lastIndex(),  k_bspline + (d+d_) - 1);
        return;
    }
    k_wavelet_last   = secondbasis.rangeJ(j_wavelet).lastIndex();
    k_wavelet_first  = k_wavelet_last - (secondbasis.cardJR(j_wavelet) + d/2) + 1;
    k_wavelet_first  = std::max((FLENS_DEFAULT_INDEXTYPE) 1, k_wavelet_first);


    return;
}

template <typename T>
void
Basis<T,Primal,Interval,Dijkema>::
getWaveletNeighborsForBSpline(FLENS_DEFAULT_INDEXTYPE j_bspline, FLENS_DEFAULT_INDEXTYPE k_bspline,
                              const Basis<T,Primal,Periodic,CDF> &secondbasis,
                              FLENS_DEFAULT_INDEXTYPE &j_wavelet, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                              FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const
{
    j_wavelet = j_bspline;

    // Translation index on intervals = periodic index + 1 // Is that right for d > 2??

    // Get wavelets for (periodic) translation index
    FLENS_DEFAULT_INDEXTYPE kminR, kmaxR;
    // Left Boundary
    if(k_bspline <= mra.rangeIL(j_bspline).lastIndex()){
    	//std::cout << "Left" << std::endl;
        k_bspline--;
        kminR = - std::ceil(0.5*(secondbasis.d + secondbasis.d_)) + 1;
        kmaxR = k_bspline + secondbasis.d - 1 + std::ceil(0.5*((secondbasis.d&1)+secondbasis.d_)) - 1;
    }
    // Right Boundary
    else if (k_bspline >= mra.rangeIR(j_bspline).firstIndex()){
    	//std::cout << "Right" << std::endl;
    	k_bspline--;
        kminR = k_bspline - secondbasis.d + std::ceil(0.5*((secondbasis.d&1)-secondbasis.d_)) + 1;
        kmaxR =  pow2i<T>(j_bspline) - 1 + std::ceil(0.5*(secondbasis.d+secondbasis.d_)) - 1;
    }
    // Inner Bspline
    else{
    	//std::cout << "Inner" << std::endl;

        k_bspline--;

        kminR = k_bspline - secondbasis.d + std::floor(0.5*((secondbasis.d&1)-secondbasis.d_)) + 1;
        kmaxR = k_bspline + secondbasis.d - 1 + std::ceil(0.5*((secondbasis.d&1)+secondbasis.d_)) - 1;
    }

    // Wrap Indizes
    k_wavelet_first = kminR >= 1? kminR : secondbasis.rangeJ(j_wavelet).lastIndex() + kminR;
    k_wavelet_last = kmaxR <= secondbasis.rangeJ(j_wavelet).lastIndex()? kmaxR : kmaxR - secondbasis.rangeJ(j_wavelet).lastIndex();

    // If we need complete basis anyway, just take consecutive index set
    if(k_wavelet_first == k_wavelet_last || k_wavelet_first == k_wavelet_last+1){
    	k_wavelet_first = secondbasis.rangeJ(j_wavelet).firstIndex();
    	k_wavelet_last = secondbasis.rangeJ(j_wavelet).lastIndex();
    }


}

template <typename T>
template <typename SecondRefinementBasis>
void
Basis<T,Primal,Interval,Dijkema>::
getBSplineNeighborsForBSpline(FLENS_DEFAULT_INDEXTYPE j_bspline1, FLENS_DEFAULT_INDEXTYPE k_bspline1,
                              const SecondRefinementBasis &secondrefinementbasis,
                              FLENS_DEFAULT_INDEXTYPE &j_bspline2, FLENS_DEFAULT_INDEXTYPE &k_bspline2_first, FLENS_DEFAULT_INDEXTYPE &k_bspline2_last) const
{
    ct_assert(SecondRefinementBasis::Side==Primal and SecondRefinementBasis::Domain==Interval
              and SecondRefinementBasis::Cons==Dijkema);
    //if (flens::IsSame<Basis<T,Primal,Interval,Dijkema>, SecondRefinementBasis>::value)

    j_bspline2 = j_bspline1;
    Support<T> supp = mra.phi.support(j_bspline1,k_bspline1);
    T a = supp.l1, b = supp.l2;
    if (a==0.L) {
        k_bspline2_first = (FLENS_DEFAULT_INDEXTYPE)secondrefinementbasis.mra.rangeI(j_bspline2).firstIndex();
        k_bspline2_last  = std::min(k_bspline1 + d, (FLENS_DEFAULT_INDEXTYPE)secondrefinementbasis.mra.rangeI(j_bspline2).lastIndex());
        return;
    }
    if (b<1.L) {
        k_bspline2_first = std::max(k_bspline1-d+1, (FLENS_DEFAULT_INDEXTYPE)secondrefinementbasis.mra.rangeI(j_bspline2).firstIndex());
        k_bspline2_last  = std::min(k_bspline1+d-1, (FLENS_DEFAULT_INDEXTYPE)secondrefinementbasis.mra.rangeI(j_bspline2).lastIndex());
        return;
    }
    k_bspline2_first = std::max((FLENS_DEFAULT_INDEXTYPE)secondrefinementbasis.mra.rangeI(j_bspline2).firstIndex(),k_bspline1 - d);
    k_bspline2_last  = (FLENS_DEFAULT_INDEXTYPE)secondrefinementbasis.mra.rangeI(j_bspline2).lastIndex();

    return;
}


template <typename T>
void
Basis<T,Primal,Interval,Dijkema>::
getScalingNeighborsForScaling(FLENS_DEFAULT_INDEXTYPE j_scaling1, FLENS_DEFAULT_INDEXTYPE k_scaling1, const Basis<T,Primal,Interval,Dijkema> &secondbasis,
                              FLENS_DEFAULT_INDEXTYPE &j_scaling2, FLENS_DEFAULT_INDEXTYPE &k_scaling_first, FLENS_DEFAULT_INDEXTYPE &k_scaling_last) const
{
    this->getBSplineNeighborsForBSpline(j_scaling1,k_scaling1,secondbasis,
                                        j_scaling2,k_scaling_first,k_scaling_last);
}

template <typename T>
void
Basis<T,Primal,Interval,Dijkema>::
getScalingNeighborsForScaling(FLENS_DEFAULT_INDEXTYPE j_scaling1, FLENS_DEFAULT_INDEXTYPE k_scaling1, const Basis<T,Primal,Periodic,CDF> &secondbasis,
                              FLENS_DEFAULT_INDEXTYPE &j_scaling2, FLENS_DEFAULT_INDEXTYPE &k_scaling_first, FLENS_DEFAULT_INDEXTYPE &k_scaling_last) const
{
    /*this->getBSplineNeighborsForBSpline(j_scaling1,k_scaling1,secondbasis,
                                        j_scaling2,k_scaling_first,k_scaling_last);
    */
    
    j_scaling2 = j_scaling1;

    // Use that k_interval = k_periodic + 1
    FLENS_DEFAULT_INDEXTYPE kminR = (k_scaling1-1) - d + 1;
    FLENS_DEFAULT_INDEXTYPE kmaxR = (k_scaling1-1) + d - 1;
    
    while(kminR < 1){
        if(secondbasis.mra.phi.phiR.support(j_scaling2, kminR).l2 > 0){
            break;
        }
        kminR++;
    }
    while(kmaxR > secondbasis.rangeJ(j_scaling2).lastIndex()){
        if(secondbasis.mra.phi.phiR.support(j_scaling2, kmaxR).l1 < 1){
            break;
        }
        kmaxR--;
    }

    k_scaling_first = kminR >= 1? kminR : (FLENS_DEFAULT_INDEXTYPE)secondbasis.mra.rangeI(j_scaling2).lastIndex() + kminR;
    k_scaling_last = kmaxR <= (FLENS_DEFAULT_INDEXTYPE)secondbasis.mra.rangeI(j_scaling2).lastIndex()? 
                                kmaxR : kmaxR - (FLENS_DEFAULT_INDEXTYPE)secondbasis.mra.rangeI(j_scaling2).lastIndex();

    if(k_scaling_first == k_scaling_last || k_scaling_first == k_scaling_last+1){
    	k_scaling_first = (FLENS_DEFAULT_INDEXTYPE)secondbasis.mra.rangeI(j_scaling2).firstIndex();
    	k_scaling_last = (FLENS_DEFAULT_INDEXTYPE)secondbasis.mra.rangeI(j_scaling2).lastIndex();
    }
}

template <typename T>
void
Basis<T,Primal,Interval,Dijkema>::
getWaveletNeighborsForScaling(FLENS_DEFAULT_INDEXTYPE j_scaling1, FLENS_DEFAULT_INDEXTYPE k_scaling1, const Basis<T,Primal,Interval,Dijkema> &secondbasis,
                              FLENS_DEFAULT_INDEXTYPE &j_wavelet, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first, FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const
{
    this->getWaveletNeighborsForBSpline(j_scaling1,k_scaling1,secondbasis,
                                        j_wavelet,k_wavelet_first,k_wavelet_last);
}

template <typename T>
void
Basis<T,Primal,Interval,Dijkema>::
getWaveletNeighborsForScaling(FLENS_DEFAULT_INDEXTYPE j_scaling1, FLENS_DEFAULT_INDEXTYPE k_scaling1, const Basis<T,Primal,Periodic,CDF> &secondbasis,
                              FLENS_DEFAULT_INDEXTYPE &j_wavelet, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first, FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const
{
    /*ct_assert(SecondBasis::Side==Primal and SecondBasis::Domain==Interval
              and SecondBasis::Cons==Dijkema);
    this->getWaveletNeighborsForBSpline(j_scaling1,k_scaling1,secondbasis,
                                        j_wavelet,k_wavelet_first,k_wavelet_last);
    */
    j_wavelet = j_scaling1;

    // Use k_interval = k_periodic + 1
    FLENS_DEFAULT_INDEXTYPE kminR = (k_scaling1-1) - d + std::floor(0.5*((d&1)-d_)) + 1;
    FLENS_DEFAULT_INDEXTYPE kmaxR = (k_scaling1-1) + d - 1 + std::ceil(0.5*((d&1)+d_)) - 1;
    
    while(kminR < 1){
        if(secondbasis.psi.psiR.support(j_wavelet, kminR).l2 > 0){
            break;
        }
        kminR++;
    }
    while(kmaxR > rangeJ(j_wavelet).lastIndex()){
        if(secondbasis.psi.psiR.support(j_wavelet, kmaxR).l1 < 1){
            break;
        }
        kmaxR--;
    }

    k_wavelet_first = kminR >= 1? kminR : secondbasis.rangeJ(j_wavelet).lastIndex() + kminR;
    k_wavelet_last = kmaxR <= secondbasis.rangeJ(j_wavelet).lastIndex()? 
                        kmaxR : kmaxR - secondbasis.rangeJ(j_wavelet).lastIndex();

    if(k_wavelet_first == k_wavelet_last || k_wavelet_first == k_wavelet_last+1){
    	k_wavelet_first = secondbasis.rangeJ(j_wavelet).firstIndex();
    	k_wavelet_last = secondbasis.rangeJ(j_wavelet).lastIndex();
    }
}

template <typename T>
template <typename SecondBasis>
void
Basis<T,Primal,Interval,Dijkema>::
getBSplineNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet, FLENS_DEFAULT_INDEXTYPE k_wavelet, const SecondBasis &secondbasis,
                              FLENS_DEFAULT_INDEXTYPE &j_bspline, FLENS_DEFAULT_INDEXTYPE &k_bspline_first, FLENS_DEFAULT_INDEXTYPE &k_bspline_last) const
{
    ct_assert(SecondBasis::Side==Primal and SecondBasis::Domain==Interval
              and SecondBasis::Cons==Dijkema);

    j_bspline = j_wavelet;
    Support<T> supp = psi.support(j_wavelet,k_wavelet);
    T a = supp.l1, b = supp.l2;
    if (a==0.) {
        k_bspline_first = secondbasis.mra.rangeI(j_bspline).firstIndex();
        k_bspline_last =  k_bspline_first+ secondbasis.mra.cardIL(j_bspline) + d + 2;
        k_bspline_last  = std::min(k_bspline_last, (FLENS_DEFAULT_INDEXTYPE)secondbasis.mra.rangeIR(j_bspline).lastIndex());
        return;
    }
    if (0.<a && b<1.) {
        k_bspline_first = std::max((FLENS_DEFAULT_INDEXTYPE)secondbasis.mra.rangeIL(j_bspline).firstIndex(), k_wavelet - (d+d_) + 1);
        k_bspline_last  = std::min((FLENS_DEFAULT_INDEXTYPE)secondbasis.mra.rangeIR(j_bspline).lastIndex(),  k_wavelet + (d+d_) - 1);
        return;
    }
    k_bspline_last   = secondbasis.mra.rangeI(j_bspline).lastIndex();
    k_bspline_first  = k_bspline_last - (secondbasis.mra.cardIR(j_bspline) + d + 2);
    k_bspline_first  = std::max((FLENS_DEFAULT_INDEXTYPE)secondbasis.mra.rangeIL(j_bspline).firstIndex(), k_bspline_first);
}

template <typename T>
template <typename SecondBasis>
void
Basis<T,Primal,Interval,Dijkema>::getScalingNeighborsForWavelet
                                  (FLENS_DEFAULT_INDEXTYPE j_wavelet, FLENS_DEFAULT_INDEXTYPE k_wavelet, const SecondBasis &secondbasis,
                                   FLENS_DEFAULT_INDEXTYPE &j_scaling, FLENS_DEFAULT_INDEXTYPE &k_scaling_first, FLENS_DEFAULT_INDEXTYPE &k_scaling_last) const
{
    ct_assert(SecondBasis::Side==Primal and SecondBasis::Domain==Interval
              and SecondBasis::Cons==Dijkema);
    this->getBSplineNeighborsForWavelet(j_wavelet,k_wavelet,secondbasis,
                                        j_scaling,k_scaling_first,k_scaling_last);

}

template <typename T>
void
Basis<T,Primal,Interval,Dijkema>::getWaveletNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet1, FLENS_DEFAULT_INDEXTYPE k_wavelet1,
                                                                const Basis<T,Primal,Interval,Dijkema> &secondbasis,
                                                                FLENS_DEFAULT_INDEXTYPE &j_wavelet2, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                                                FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const
{

    j_wavelet2 = j_wavelet1;
    Support<T> supp = psi.support(j_wavelet1,k_wavelet1);
    T a = supp.l1, b = supp.l2;
    if (a==0.L) {
        k_wavelet_first = (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).firstIndex();
        k_wavelet_last  = std::min(k_wavelet1 + 2*d, (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).lastIndex());
        return;
    }
    if (b<1.L) {
        k_wavelet_first = std::max(k_wavelet1-2*d+1, (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).firstIndex());
        k_wavelet_last  = std::min(k_wavelet1+2*d-1, (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).lastIndex());
        return;
    }
    k_wavelet_first = std::max((FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).firstIndex(),k_wavelet1 - 2*d);
    k_wavelet_last  = (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).lastIndex();

    return;
}

template <typename T>
void
Basis<T,Primal,Interval,Dijkema>::getWaveletNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet1, FLENS_DEFAULT_INDEXTYPE k_wavelet1,
                                                                const Basis<T,Primal,Periodic,CDF> &secondbasis,
                                                                FLENS_DEFAULT_INDEXTYPE &j_wavelet2, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                                                FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const
{
    j_wavelet2 = j_wavelet1;

    // Use that k_interval = k_periodic + 1
    FLENS_DEFAULT_INDEXTYPE kminR = (k_wavelet1-1) - d - d_ + 1 + 1;
    FLENS_DEFAULT_INDEXTYPE kmaxR = (k_wavelet1-1) + d + d_ - 1 - 1;
    
    while(kminR < 1){
        if(secondbasis.psi.psiR.support(j_wavelet2, kminR).l2 > 0){
            break;
        }
        kminR++;
    }
    while(kmaxR > rangeJ(j_wavelet2).lastIndex()){
        if(secondbasis.psi.psiR.support(j_wavelet2, kmaxR).l1 < 1){
            break;
        }
        kmaxR--;
    }

    k_wavelet_first = kminR >= 1? kminR : secondbasis.rangeJ(j_wavelet2).lastIndex() + kminR;
    k_wavelet_last = kmaxR <= secondbasis.rangeJ(j_wavelet2).lastIndex()? kmaxR : kmaxR - secondbasis.rangeJ(j_wavelet2).lastIndex();

    if(k_wavelet_first == k_wavelet_last || k_wavelet_first == k_wavelet_last+1){
    	k_wavelet_first = secondbasis.rangeJ(j_wavelet2).firstIndex();
    	k_wavelet_last = secondbasis.rangeJ(j_wavelet2).lastIndex();
    }

}

template <typename T>
template <typename SecondBasis>
void
Basis<T,Primal,Interval,Dijkema>::getLowerWaveletNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet1, FLENS_DEFAULT_INDEXTYPE k_wavelet1,
                                                                     const SecondBasis &secondbasis,
                                                                     FLENS_DEFAULT_INDEXTYPE &j_wavelet2, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                                                     FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const
{
    ct_assert(SecondBasis::Side==Primal and SecondBasis::Domain==Interval
              and SecondBasis::Cons==Dijkema);
    //if (flens::IsSame<Basis<T,Primal,Interval,Dijkema>, SecondRefinementBasis>::value)

    j_wavelet2 = j_wavelet1-1;
    Support<T> supp = psi.support(j_wavelet1,k_wavelet1);
    T a = supp.l1, b = supp.l2;
    FLENS_DEFAULT_INDEXTYPE k_tilde = k_wavelet1/2;
    if (a==0.L) {
        k_wavelet_first = (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).firstIndex();
        k_wavelet_last  = std::min(k_tilde + 2*d, (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).lastIndex());
        return;
    }
    if (b<1.L) {
        k_wavelet_first = std::max(k_tilde-2*d+1, (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).firstIndex());
        k_wavelet_last  = std::min(k_tilde+2*d-1, (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).lastIndex());
        return;
    }
    k_wavelet_first = std::max((FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).firstIndex(),k_tilde - 2*d);
    k_wavelet_last  = (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).lastIndex();

    return;
}

template <typename T>
void
Basis<T,Primal,Interval,Dijkema>::getHigherWaveletNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet1, FLENS_DEFAULT_INDEXTYPE k_wavelet1,
                                                                     const Basis<T,Primal,Interval,Dijkema> &secondbasis,
                                                                     FLENS_DEFAULT_INDEXTYPE &j_wavelet2, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                                                     FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const
{

    j_wavelet2 = j_wavelet1+1;
    Support<T> supp = psi.support(j_wavelet1,k_wavelet1);
    T a = supp.l1, b = supp.l2;
    FLENS_DEFAULT_INDEXTYPE k_tilde = k_wavelet1*2;
    if (a==0.L) {
        k_wavelet_first = (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).firstIndex();
        k_wavelet_last  = std::min(k_tilde + 3*d, (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).lastIndex());
        return;
    }
    if (b<1.L) {
        k_wavelet_first = std::max(k_tilde-2*d-1, (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).firstIndex());
        k_wavelet_last  = std::min(k_tilde+2*d+1, (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).lastIndex());
        return;
    }
    k_wavelet_first = std::max((FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).firstIndex(),k_tilde - 3*d);
    k_wavelet_last  = (FLENS_DEFAULT_INDEXTYPE)secondbasis.rangeJ(j_wavelet2).lastIndex();

    return;
}

template <typename T>
void
Basis<T,Primal,Interval,Dijkema>::getHigherWaveletNeighborsForWavelet(FLENS_DEFAULT_INDEXTYPE j_wavelet1, FLENS_DEFAULT_INDEXTYPE k_wavelet1,
                                                                     const Basis<T,Primal,Periodic,CDF> &secondbasis,
                                                                     FLENS_DEFAULT_INDEXTYPE &j_wavelet2, FLENS_DEFAULT_INDEXTYPE &k_wavelet_first,
                                                                     FLENS_DEFAULT_INDEXTYPE &k_wavelet_last) const
{    
    j_wavelet2 = j_wavelet1+1;
    // Use that k_interval = k_periodic + 1

    FLENS_DEFAULT_INDEXTYPE kminR = 2 + 2*(k_wavelet1-1) - std::ceil(1.5*(d + d_)) + 1;
    FLENS_DEFAULT_INDEXTYPE kmaxR = 2*(k_wavelet1-1) - 1 + std::ceil(1.5*(d+d_)) - 1;

    while(kminR < 1){
        if(secondbasis.psi.psiR.support(j_wavelet2, kminR).l2 > 0){
            break;
        }
        kminR++;
    }
    while(kmaxR > rangeJ(j_wavelet2).lastIndex()){
        if(secondbasis.psi.psiR.support(j_wavelet2, kmaxR).l1 < 1){
            break;
        }
        kmaxR--;
    }
    
    k_wavelet_first = kminR >= 1? kminR : secondbasis.rangeJ(j_wavelet2).lastIndex() + kminR;
    k_wavelet_last = kmaxR <= secondbasis.rangeJ(j_wavelet2).lastIndex()? 
                    kmaxR : kmaxR - secondbasis.rangeJ(j_wavelet2).lastIndex();

    if(k_wavelet_first == k_wavelet_last || k_wavelet_first == k_wavelet_last+1){
    	k_wavelet_first = secondbasis.rangeJ(j_wavelet2).firstIndex();
    	k_wavelet_last = secondbasis.rangeJ(j_wavelet2).lastIndex();
    }
    
    return;
}

template <typename T>
Basis<T,Primal,Interval,Dijkema>::LaplaceOperator1D::
LaplaceOperator1D(FLENS_DEFAULT_INDEXTYPE _d, const Basis<T,Primal,Interval,Dijkema> &_refinementbasis)
 : d(_d), refinementbasis(_refinementbasis)
{
    switch (d) {
        case 2:
            inner_values.engine().resize((FLENS_DEFAULT_INDEXTYPE)2,0);
            inner_values = 2., -1.;
            break;
        case 3:
            outer_values.engine().resize((FLENS_DEFAULT_INDEXTYPE)3,0);
            outer_values = 4.L/3.L, -1.L/6.L, -1.L/6.L;
            inner_values.engine().resize((FLENS_DEFAULT_INDEXTYPE)3,0);
            inner_values = 1.L, -1.L/3.L, -1.L/6.L;
            break;
        default:
            std::cerr << "Basis<T,Primal,Interval,Dijkema>::LaplaceOperator1D "
                      << "not yet implemented for d=" << d << std::endl;
    }
}

template <typename T>
T
Basis<T,Primal,Interval,Dijkema>::LaplaceOperator1D::
operator()(XType xtype1, FLENS_DEFAULT_INDEXTYPE j1, FLENS_DEFAULT_INDEXTYPE k1, XType xtype2, FLENS_DEFAULT_INDEXTYPE j2, FLENS_DEFAULT_INDEXTYPE k2)
{
    assert(j1==j2);
    assert(refinementbasis._bc(0) == 1);
    assert(refinementbasis._bc(1) == 1);

    if (d==2) {
        FLENS_DEFAULT_INDEXTYPE k_diff = std::abs(k1 - k2);
        if (k_diff>1) return 0.L;
        return  pow2i<T>(2*j1)*inner_values(k_diff);
    }
    else if (d==3) {
        FLENS_DEFAULT_INDEXTYPE k_diff = std::abs(k1 - k2);
        if (k_diff>2) return 0.L;

        FLENS_DEFAULT_INDEXTYPE firstIndex =  refinementbasis.mra.rangeI(j1).firstIndex();
        FLENS_DEFAULT_INDEXTYPE lastIndex  =  refinementbasis.mra.rangeI(j1).lastIndex();
        if (k1==firstIndex || k1==lastIndex || k2 == firstIndex || k2 == lastIndex) {
            return pow2i<T>(2*j1)*outer_values(k_diff);
        }
        else {
            return pow2i<T>(2*j1)*inner_values(k_diff);
        }
    }
    else {
        std::cerr << "Basis<T,Primal,Interval,Dijkema>::PoissonOperator1D::operator() "
                  << "not yet implemented for d=" << d << std::endl;
        exit(1);
        return 0.;
    }
}

template <typename T>
Basis<T,Primal,Interval,Dijkema>::IdentityOperator1D::
IdentityOperator1D(FLENS_DEFAULT_INDEXTYPE _d, const Basis<T,Primal,Interval,Dijkema> &_refinementbasis)
 : d(_d), refinementbasis(_refinementbasis)
{
    switch (d) {
        case 2:
            inner_values.engine().resize((FLENS_DEFAULT_INDEXTYPE)2,0);
            inner_values = 2.L/3.L, 1.L/6.L;
            break;
        case 3:
            outer_values.engine().resize((FLENS_DEFAULT_INDEXTYPE)3,0);
            outer_values = 1.L/3.L, 5.L/24.L, 1.L/120.L;
            inner_values.engine().resize((FLENS_DEFAULT_INDEXTYPE)3,0);
            inner_values = 11.L/20.L, 13.L/60.L, 1.L/120.L;
            break;
        default:
            std::cerr << "Basis<T,Primal,Interval,Dijkema>::IdentityOperator1D "
                      << "not yet implemented for d=" << d << std::endl;
    }
}

template <typename T>
T
Basis<T,Primal,Interval,Dijkema>::IdentityOperator1D::
operator()(XType xtype1, FLENS_DEFAULT_INDEXTYPE j1, FLENS_DEFAULT_INDEXTYPE k1, XType xtype2, FLENS_DEFAULT_INDEXTYPE j2, FLENS_DEFAULT_INDEXTYPE k2)
{
    assert(j1==j2);
    assert(refinementbasis._bc(0) == 1);
    assert(refinementbasis._bc(1) == 1);
    if (d==2) {
        FLENS_DEFAULT_INDEXTYPE k_diff = std::abs(k1 - k2);
        if (k_diff>1) return 0.L;
        return  inner_values(k_diff);
    }
    else if (d==3) {
        FLENS_DEFAULT_INDEXTYPE k_diff = std::abs(k1 - k2);
        if (k_diff>2) return 0.L;

        FLENS_DEFAULT_INDEXTYPE firstIndex =  refinementbasis.mra.rangeI(j1).firstIndex();
        FLENS_DEFAULT_INDEXTYPE lastIndex  =  refinementbasis.mra.rangeI(j1).lastIndex();
        if (k1==firstIndex || k1==lastIndex || k2 == firstIndex || k2 == lastIndex) {
            return outer_values(k_diff);
        }
        else {
            return inner_values(k_diff);
        }
    }
    else {
        std::cerr << "Basis<T,Primal,Interval,Dijkema>::IdentityOperator1D::operator() "
                  << "not yet implemented for d=" << d << std::endl;
        exit(1);
        return 0.;
    }
}

} // namespace lawa

