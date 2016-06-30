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

#ifndef LAWA_INTEGRALS_EXPWEIGHTINTEGRAL_H
#define LAWA_INTEGRALS_EXPWEIGHTINTEGRAL_H 1

#include <lawa/settings/enum.h>
#include <lawa/functiontypes/exponentialweightfunction1d.h>
#include <lawa/integrals/quadrature.h>

namespace lawa {

template <QuadratureType Quad, typename First, typename Second=First>
struct IntegralExpWeight
{
    typedef typename First::T T;

    IntegralExpWeight(const First &_first, const T eta=2., const T _left=0., const T _right=1.);

    IntegralExpWeight(const First &_first, const Second &_second, const T eta=2.,
                      const T _left=0., const T _right=1.);

    T
    operator()(FLENS_DEFAULT_INDEXTYPE _j1, FLENS_DEFAULT_INDEXTYPE _k1, XType _e1, FLENS_DEFAULT_INDEXTYPE _deriv) const;

    T
    operator()(FLENS_DEFAULT_INDEXTYPE _j1, FLENS_DEFAULT_INDEXTYPE _k1, XType _e1, FLENS_DEFAULT_INDEXTYPE _deriv1,
               FLENS_DEFAULT_INDEXTYPE _j2, FLENS_DEFAULT_INDEXTYPE _k2, XType _e2, FLENS_DEFAULT_INDEXTYPE _deriv2) const;

    T
    integrand(T x) const;

    ExponentialWeightFunction1D<T> expweight;
    const First &first;
    const Second &second;
    const T left;
    const T right;
    const T RightmLeft;
    const T SqrtRightmLeft;
    const T OneDivSqrtRightmLeft;
    Quadrature<Quad,IntegralExpWeight<Quad,First,Second> > quadrature;
    const bool _f2;
    mutable FLENS_DEFAULT_INDEXTYPE j1, deriv1,
                j2, deriv2;
    mutable FLENS_DEFAULT_INDEXTYPE k1, k2;
    mutable XType e1, e2;
};

} // namespace lawa

#include <lawa/integrals/expweightintegral.tcc>

#endif // LAWA_INTEGRALS_EXPWEIGHTINTEGRAL_H

