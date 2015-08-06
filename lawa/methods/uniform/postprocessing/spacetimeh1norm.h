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
 
#ifndef LAWA_METHODS_UNIFORM_POSTPROCESSING_SPACETIMEH1NORM_H
#define LAWA_METHODS_UNIFORM_POSTPROCESSING_SPACETIMEH1NORM_H 1

namespace lawa {

template <typename T, typename Basis2D>
class SpaceTimeH1Norm
{
    private:
        typedef typename Basis2D::FirstBasisType Basis_t;
        typedef typename Basis2D::SecondBasisType Basis_x;

        Integral<Gauss, Basis_t, Basis_t> integral_t;
        Integral<Gauss, Basis_x, Basis_x> integral_x; 

    public:

        const Basis2D &basis;

        SpaceTimeH1Norm(const Basis2D &_basis);

        /* Calculates the H1(0,T; H1)-norm (or W(0,T)-norm) for u = u1(t)*u2(x),
         *  ||u|| =  sqrt( ||u1||_L2 * ||u2||_H1 + |u1|_H1 * ||u2||_(H1)' )
         */
        T
        operator()(XType xtype_t, int j_t, int k_t, 
                   XType xtype_x, int j_x, int k_x) const;
};

} // namespace lawa

#include <lawa/methods/uniform/postprocessing/spacetimeh1norm.tcc>

#endif // LAWA_METHODS_UNIFORM_POSTPROCESSING_SPACETIMEH1NORM_H

