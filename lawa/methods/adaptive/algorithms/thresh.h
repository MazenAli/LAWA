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
 
#ifndef  LAWA_METHODS_ADAPTIVE_ALGORITHMS_THRESH_H
#define  LAWA_METHODS_ADAPTIVE_ALGORITHMS_THRESH_H 1

#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/datastructures/index.h>

namespace lawa {

template <typename T, typename Index>
    Coefficients<Lexicographical,T,Index >
    ABSOLUTE_THRESH(const Coefficients<Lexicographical,T,Index > &v, T eta);

template <typename T>
    Coefficients<Lexicographical,T,Index1D >
    THRESH(const Coefficients<Lexicographical,T,Index1D > &v, T eta, bool deleteBSpline=true,
           bool hp=false);

template <typename T>
    Coefficients<Lexicographical,T,Index2D >
    THRESH(const Coefficients<Lexicographical,T,Index2D > &v, T eta, bool deleteBSpline=true,
           bool hp=false);

template <typename T>
    void
    THRESH_NoCopy(Coefficients<Lexicographical,T,Index2D > &v, T eta);

template <typename T>
    Coefficients<Lexicographical,T,Index2D >
    MULTITREE_THRESH(const Coefficients<Lexicographical,T,Index2D > &v, T eta);

template <typename T>
    Coefficients<Lexicographical,T,Index3D >
    THRESH(const Coefficients<Lexicographical,T,Index3D > &v, T eta, bool deleteBSpline=true,
           bool hp=false);

template <typename T, typename Index>
    Coefficients<Lexicographical,T,Index >
    THRESH(const Coefficients<AbsoluteValue,T,Index > &v, T eta, bool hp=false);



} // namespace lawa

#include <lawa/methods/adaptive/algorithms/thresh.tcc>

#endif // LAWA_METHODS_ADAPTIVE_ALGORITHMS_THRESH_H

