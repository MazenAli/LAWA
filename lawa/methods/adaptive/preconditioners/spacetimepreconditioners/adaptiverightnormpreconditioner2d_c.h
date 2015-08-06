/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2014  Sebastian Kestler, Mario Rometsch, Kristina Steih, 
  Alexander Stippler.

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

#ifndef LAWA_METHODS_ADAPTIVE_PRECONDITIONERS_SPACETIMEPRECONDTIONERS_ADAPTIVERIGHTNORMPRECONDITIONER2D_H
#define LAWA_METHODS_ADAPTIVE_PRECONDITIONERS_SPACETIMEPRECONDTIONERS_ADAPTIVERIGHTNORMPRECONDITIONER2D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/integrals/integrals.h>
#include <lawa/settings/enum.h>

namespace lawa {

/* Right Preconditioner for Space-Time Problems
 *     Normalization in  X = L_2(0,T) \otimes V   \cap   H_1(0,t) \otimes V' (e.g. V = H^1):
 *          Prec = ( || bf_t ||_L2^2 * || bf_x ||_V^2 + || bf_t ||_H1^2 * || bf_x ||_V'^2 )^(-1/2)
 *     (see Schwab/Stevenson "Adaptive Wavelet Methods for Parabolic Problems", Math.Comp. 78 (267), pp. 1293-1318, 2009)
 *
 *	Computations:
 *		|| bf_t ||_L2^2 : Integral								   ( s!=2 -> 1)
 *		|| bf_t ||_H1^2 : Scaling of Integral
 *
 *		|| bf_x ||_H1^2 : Integral 								   ( s!=2 -> Scaling of Integral)
 *		|| bf_x ||_H-1^2: Duality ( ||.||_H-1 \equiv 1./||.||_H1)  ( s!=2 -> Scaling of 1)
 *
 *  Adaptive: Computes values are stored in a hashmap
 */
template <typename T, typename Basis2D>
class AdaptiveRightNormPreconditioner2D_c
{
    typedef typename Basis2D::FirstBasisType  FirstBasis;
    typedef typename Basis2D::SecondBasisType SecondBasis;
    typedef Coefficients<Lexicographical,T,Index1D>		Coeffs_1D;

    public:
        AdaptiveRightNormPreconditioner2D_c(const Basis2D &basis, T s=2.);  //s=2: A: H^1 -> H^{-1}

        T
        operator()(XType xtype1, int j1, long k1,
                   XType xtype2, int j2, long k2);

        T
        operator()(const Index2D &index);

    private:
        //const Basis2D &_basis;
        T              _s;  // scaling for certain classes of integral operators
        Integral<Gauss,FirstBasis,FirstBasis>   _integral_t;
        Integral<Gauss,SecondBasis,SecondBasis> _integral_x;

        Coeffs_1D values_L2_t, values_L2_x, values_H1semi;
};

}   // namespace lawa

#include <lawa/methods/adaptive/preconditioners/spacetimepreconditioners/adaptiverightnormpreconditioner2d_c.tcc>

#endif // LAWA_METHODS_ADAPTIVE_PRECONDITIONERS_SPACETIMEPRECONDTIONERS_ADAPTIVERIGHTNORMPRECONDITIONER2D_H
