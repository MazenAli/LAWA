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

#ifndef LAWA_FUNCTIONTYPES_EXPONENTIALWEIGHTFUNCTION1D_H
#define LAWA_FUNCTIONTYPES_EXPONENTIALWEIGHTFUNCTION1D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/settings/enum.h>

namespace lawa {

template<typename T>
struct ExponentialWeightFunction1D
{
    static T eta;
    static T x1, x2;
    static bool with_precPoints;
    static flens::DenseVector<flens::Array<T> > singularPoints;

    static void
    setEta(T _eta);

    static void
    setSingularPoints();

    static void
    setPrecPoints(T _x1, T _x2);

    static T
    alpha(T x);

    static T
    dalpha(T x);


    static T
    weight(T x);

    static T
    sqrtweight(T x);

    static T
    dweight(T x);
};


}   //namespace lawa

#include <lawa/functiontypes/exponentialweightfunction1d.tcc>

#endif  // LAWA_FUNCTIONTYPES_EXPONENTIALWEIGHTFUNCTION1D_H
