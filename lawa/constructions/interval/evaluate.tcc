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

#include <lawa/aux/compiletime_assert.h>

namespace lawa {

template <FunctionSide Side, Construction Cons, typename X>
typename X::ElementType
evaluate(const MRA<typename X::ElementType,Side,Interval,Cons> &mra, FLENS_DEFAULT_INDEXTYPE j,
         const flens::DenseVector<X>& coeffs, typename X::ElementType x, FLENS_DEFAULT_INDEXTYPE deriv)
{
    ct_assert(Side==Primal or Side==Orthogonal);

    assert(j>=mra.j0);
    assert(coeffs.length()==mra.cardI(j));
    assert(x>=0.);
    assert(x<=1.);

    typedef typename X::ElementType T;

    BSpline<T,Side,Interval,Cons> phi(mra);
    T ret = 0.0;
    for (FLENS_DEFAULT_INDEXTYPE k=mra.rangeI(j).firstIndex(); k<=mra.rangeI(j).lastIndex(); ++k) {
        ret += coeffs(k) * phi(x,j,k,deriv);
    }
    return ret;
/*    if (x<mra.suppIL(j).l2) {
        for (FLENS_DEFAULT_INDEXTYPE k=mra.rangeIL(j).firstIndex(); k<=mra.rangeIL(j).lastIndex(); ++k) {
            ret += coeffs(k) * phi(x,j,k);
        }
    }
    if (x>mra.suppIR(j).l1) {
        for (FLENS_DEFAULT_INDEXTYPE k=mra.rangeIR(j).firstIndex(); k<=mra.rangeIR(j).lastIndex(); ++k) {
            ret += coeffs(k) * phi(x,j,k);
        }
    }

    FLENS_DEFAULT_INDEXTYPE shift = ifloor(pow2i<T>(j) * fabs(mra.phi.support(j,mra.rangeII(j).firstIndex()).l1-x));
    FLENS_DEFAULT_INDEXTYPE firstk = mra.rangeII(j).firstIndex();
    FLENS_DEFAULT_INDEXTYPE sectork = firstk + shift;

    FLENS_DEFAULT_INDEXTYPE start = std::max(sectork-mra.d+1, firstk);
    FLENS_DEFAULT_INDEXTYPE end = std::min(sectork+mra.d-1, mra.rangeII(j).lastIndex());
    for (FLENS_DEFAULT_INDEXTYPE k=start; k<=end; ++k) {
        ret += coeffs(k) * phi(x,j,k);
    }
    return ret;
*/
}

template <FunctionSide Side, Construction Cons, typename X>
typename X::ElementType
evaluate(const Basis<typename X::ElementType,Side,Interval,Cons> &basis,
         FLENS_DEFAULT_INDEXTYPE J, const flens::DenseVector<X> &coeffs, typename X::ElementType x,
         FLENS_DEFAULT_INDEXTYPE deriv)
{
    ct_assert(Side==Primal or Side==Orthogonal);

    assert(J>=basis.j0);
    assert(coeffs.range()==basis.mra.rangeI(J));
    assert(x>=0.);
    assert(x<=1.);

    typedef typename X::ElementType T;
/*
    const FLENS_DEFAULT_INDEXTYPE j0 = basis.j0;

    basis.setLevel(j0);
    T ret = 0.;
    ret += evaluate(basis.mra,j0,coeffs(basis.mra.rangeI(j0)),x,deriv);
    Wavelet<T,Side,Interval,Cons> psi(basis);
    for (FLENS_DEFAULT_INDEXTYPE j=j0; j<=J-1; ++j) {
        if (x<basis.suppJL(j).l2) {
            for (FLENS_DEFAULT_INDEXTYPE k=basis.rangeJL(j).firstIndex(); k<=basis.rangeJL(j).lastIndex(); ++k) {
                ret += coeffs(basis.rangeI(j).lastIndex() + k) * psi(x, j, k, deriv);
            }
        }
        if (x>basis.suppJR(j).l1) {
            for (FLENS_DEFAULT_INDEXTYPE k=basis.rangeJR(j).firstIndex(); k<=basis.rangeJR(j).lastIndex(); ++k) {
                ret += coeffs(basis.rangeI(j).lastIndex() + k) * psi(x, j, k);
            }
        }

        T first = basis.suppJI(j).l1;
        T  last = basis.suppJI(j).l2;

        if ((x>first) && (x<last)) {
            FLENS_DEFAULT_INDEXTYPE shift = ifloor(pow2i<T>(j) * (x-first));
            FLENS_DEFAULT_INDEXTYPE pos = basis.rangeJI(j).firstIndex() + shift;

            FLENS_DEFAULT_INDEXTYPE from = std::max(basis.rangeJI(j).firstIndex(),pos-(basis.d+basis.d_));
            FLENS_DEFAULT_INDEXTYPE to   = std::min(basis.rangeJI(j).lastIndex(), pos+(basis.d+basis.d_));
            for (FLENS_DEFAULT_INDEXTYPE k=from; k<=to; ++k) {
                ret += coeffs(basis.rangeI(j).lastIndex()+k) * psi(x,j,k);
            }
        }
    }
    return ret;
*/
    const FLENS_DEFAULT_INDEXTYPE j0 = basis.j0;
    basis.setLevel(j0);

    T ret = 0;
    ret += evaluate(basis.mra,j0,coeffs(flens::Range<FLENS_DEFAULT_INDEXTYPE>((FLENS_DEFAULT_INDEXTYPE)basis.mra.rangeI(j0).firstIndex(),(FLENS_DEFAULT_INDEXTYPE)basis.mra.rangeI(j0).lastIndex())),x,deriv);
    Wavelet<T,Side,Interval,Cons> psi(basis);
    for (FLENS_DEFAULT_INDEXTYPE j=j0; j<=J-1; ++j) {
        for (FLENS_DEFAULT_INDEXTYPE k=1; k<=basis.cardJ(j); ++k) {
            ret += coeffs(basis.mra.rangeI(j).lastIndex() + k) * psi(x, j, k, deriv);
        }
    }
    return ret;
}

} // namespace lawa

