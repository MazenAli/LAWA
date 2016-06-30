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

template <typename T, Construction Cons>
Wavelet<T,Primal,Interval,Cons>::Wavelet(const Basis<T,Primal,Interval,Cons> &_basis)
    : basis(_basis)
{
}

template <typename T, Construction Cons>
T
Wavelet<T,Primal,Interval,Cons>::operator()(T x, FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k, unsigned short deriv) const
{
    using flens::_;

    assert(x>=0.);
    assert(x<=1.);
    assert(j>=basis.min_j0);
    assert(k>=basis.rangeJ(j).firstIndex());
    assert(k<=basis.rangeJ(j).lastIndex());

    const typename flens::DenseVector<flens::Array<T> >::ConstView coeffs = basis.M1(j,_,k);
    T ret = 0;
    for (FLENS_DEFAULT_INDEXTYPE r=coeffs.firstIndex(); r<=coeffs.lastIndex(); ++r) {
        ret += coeffs(r) * basis.mra.phi(x,j+1,r,deriv);
    }
    return ret;
}

template <typename T, Construction Cons>
Support<T>
Wavelet<T,Primal,Interval,Cons>::support(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    assert(j>=basis.min_j0);
    assert(k>=basis.rangeJ(j).firstIndex());
    assert(k<=basis.rangeJ(j).lastIndex());
    if (k<=basis.M1.left.lastIndex()) {
        return Support<T>(0.,pow2i<T>(-j-1)*basis.M1.lengths(k));
    }
    if (k>pow2i<T>(j)-basis.M1.right.length()) {
        return Support<T>(1-pow2i<T>(-j-1)
                        *(basis.M1.lengths(k-1-pow2i<T>(j))), 1.);
    }
    // FIXME: remove std::max: left support end cannot be less than 0. Check for error (Primbs!!!)
    return pow2i<T>(-j-1)*Support<T>(std::max((FLENS_DEFAULT_INDEXTYPE) 0,basis.M1.lengths(0)+1-basis.d+2*(k-basis.M1.left.lastIndex()-1)),
                                     basis.M1.lengths(0)+basis.M1.leftband.length()+2*(k-basis.M1.left.lastIndex()-1));
}

template <typename T, Construction Cons>
flens::DenseVector<flens::Array<T> >
Wavelet<T,Primal,Interval,Cons>::singularSupport(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    assert(j>=basis.min_j0);
    assert(k>=basis.rangeJ(j).firstIndex());
    assert(k<=basis.rangeJ(j).lastIndex());

    if (k<=basis.M1.left.lastIndex()) {
        return linspace((T)0.,
                        pow2i<T>(-j-1)*basis.M1.lengths(k),
                        2*basis.M1.lengths(k)+1.);
    }
    if (k>pow2i<T>(j)-basis.M1.right.length()) {
        return linspace(1-pow2i<T>(-j-1)*(basis.M1.lengths(k-1-pow2i<T>(j))),
                        (T)1.,
                        2*basis.M1.lengths(k-1-pow2i<T>(j))+1.);
    }
    // FIXME: remove std::max: left support end cannot be less than 0. Check for error (Primbs!!!)
    return pow2i<T>(-j-1)*linspace(std::max(0.,basis.M1.lengths(0)+1-basis.d+2*(k-basis.M1.left.lastIndex()-1.)),
                                   basis.M1.lengths(0)+basis.M1.leftband.length()+2*(k-basis.M1.left.lastIndex()-1.),
                                   // FIXME: understand why value in last line is too large
                                   2*(basis.d+basis.d_)-1);
                                   // FIXME: understand why 2*n instead of 2*n-1  ... (+d+1)
                                   //2*(basis.M1.leftband.length())+basis.d+1.);
}

template <typename T, Construction Cons>
FLENS_DEFAULT_INDEXTYPE
Wavelet<T,Primal,Interval,Cons>::vanishingMoments(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    assert(j>=basis.min_j0);
    assert(k>=basis.rangeJ(j).firstIndex());
    assert(k<=basis.rangeJ(j).lastIndex());

    assert(0);
    return 0;
}

template <typename T, Construction Cons>
T
Wavelet<T,Primal,Interval,Cons>::tic(FLENS_DEFAULT_INDEXTYPE j) const
{
    return pow2i<T>(-j-1);
}

template <typename T, Construction Cons>
flens::DenseVector<flens::Array<long double> > *
Wavelet<T,Primal,Interval,Cons>::getRefinement(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k, FLENS_DEFAULT_INDEXTYPE &refinement_j, FLENS_DEFAULT_INDEXTYPE &refinement_k_first,
												FLENS_DEFAULT_INDEXTYPE &split, FLENS_DEFAULT_INDEXTYPE &refinement_k_restart) const
{
	// No split necessary, so set default values
	refinement_k_restart = 1;

	k -= 1;
    refinement_j = j + 1;
    // left boundary
    if (k<basis.cardJL(j)) {
        refinement_k_first = basis._leftOffsets[k];
    	split = basis._leftRefCoeffs[k].length() + 1;
        return &(basis._leftRefCoeffs[k]);
    }
    // inner part
    if (k<basis.cardJL(j)+basis.cardJI(j)) {
        FLENS_DEFAULT_INDEXTYPE type = 0;
        if (basis.d % 2 != 0 && k+1>basis.cardJ(j)/2.) type = 1;
        FLENS_DEFAULT_INDEXTYPE shift = k+(FLENS_DEFAULT_INDEXTYPE) 1;
        refinement_k_first = 2*shift+basis._innerOffsets[type];
        return &(basis._innerRefCoeffs[type]);
    }
    // right part
    FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)(k+1 - (basis.cardJ(j) - basis.cardJL(j) + 1));
    FLENS_DEFAULT_INDEXTYPE shift = pow2i<FLENS_DEFAULT_INDEXTYPE>(j)-1;
    refinement_k_first = 2*shift+basis._rightOffsets[type];
    return &(basis._rightRefCoeffs[type]);
}

template <typename T, Construction Cons>
FLENS_DEFAULT_INDEXTYPE
Wavelet<T,Primal,Interval,Cons>::getRefinementLevel(FLENS_DEFAULT_INDEXTYPE j) const
{
    return j + 1;
}

template <typename T, Construction Cons>
void
Wavelet<T,Primal,Interval,Cons>::getRefinementNeighbors(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k, FLENS_DEFAULT_INDEXTYPE &refinement_j,
                                                        FLENS_DEFAULT_INDEXTYPE &refinement_k_first,
                                                        FLENS_DEFAULT_INDEXTYPE &refinement_k_last) const
{
    k -= 1;
    refinement_j = j + 1;
    // left boundary
    if (k<basis.cardJL(j)) {
        refinement_k_first = basis._leftOffsets[k];
        refinement_k_last  = refinement_k_first + basis._leftRefCoeffs[k].lastIndex();
        return;
    }
    // inner part
    if (k<basis.cardJL(j)+basis.cardJI(j)) {
        FLENS_DEFAULT_INDEXTYPE type = 0;
        if (basis.d % 2 != 0 && k+1>basis.cardJ(j)/2.) type = 1;
        FLENS_DEFAULT_INDEXTYPE shift = k+(FLENS_DEFAULT_INDEXTYPE) 1;
        refinement_k_first = 2*shift+basis._innerOffsets[type];
        refinement_k_last  = refinement_k_first + basis._innerRefCoeffs[type].lastIndex();
        return;
    }
    // right part
    FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)(k+1 - (basis.cardJ(j) - basis.cardJL(j) + 1));
    FLENS_DEFAULT_INDEXTYPE shift = pow2i<FLENS_DEFAULT_INDEXTYPE>(j)-1;
    refinement_k_first = 2*shift+basis._rightOffsets[type];
    refinement_k_last  = refinement_k_first + basis._rightRefCoeffs[type].lastIndex();
    return;
}

template <typename T, Construction Cons>
T
Wavelet<T,Primal,Interval,Cons>::getL2Norm(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    k -= 1;
    // left boundary
    if (k<basis.cardJL(j)) {
        return basis._leftL2Norms[k];
    }
    // inner part
    if (k<basis.cardJL(j)+basis.cardJI(j)) {
        return basis._innerL2Norms[0];
    }
    // right part
    FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)(k+1 - (basis.cardJ(j) - basis.cardJL(j) + 1));
    return basis._rightL2Norms[type];
}

template <typename T, Construction Cons>
T
Wavelet<T,Primal,Interval,Cons>::getH1SemiNorm(FLENS_DEFAULT_INDEXTYPE j, FLENS_DEFAULT_INDEXTYPE k) const
{
    T pow2ij = (T)((FLENS_DEFAULT_INDEXTYPE) 1 << j);
    k -= 1;
    // left boundary
    if (k<basis.cardJL(j)) {
        return pow2ij*basis._leftH1SemiNorms[k];
    }
    // inner part
    if (k<basis.cardJL(j)+basis.cardJI(j)) {
        return pow2ij*basis._innerH1SemiNorms[0];
    }
    // right part
    FLENS_DEFAULT_INDEXTYPE type  = (FLENS_DEFAULT_INDEXTYPE)(k+1 - (basis.cardJ(j) - basis.cardJL(j) + 1));
    return pow2ij*basis._rightH1SemiNorms[type];
}

} // namespace lawa

