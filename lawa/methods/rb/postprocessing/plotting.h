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
 
#ifndef LAWA_METHODS_RB_POSTPROCESSING_PLOTTING_H
#define LAWA_METHODS_RB_POSTPROCESSING_PLOTTING_H

#include <lawa/methods/adaptive/datastructures/datastructures.h>

namespace lawa {

// Version NonPeriodic x NonPeriodic
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<!IsPeriodic<typename Basis::FirstBasisType>::value
					and !IsPeriodic<typename Basis::SecondBasisType>::value, T>::value, void>::Type
saveCoeffVector2D(const Coefficients<Lexicographical,T,Index2D> &coeff, const Basis &basis2d, const char* filename);

// Version Periodic x NonPeriodic
template <typename T, typename Basis>
typename cxxblas::RestrictTo<SFINAE_Wrapper<IsPeriodic<typename Basis::FirstBasisType>::value
					and !IsPeriodic<typename Basis::SecondBasisType>::value, T>::value, void>::Type
saveCoeffVector2D(const Coefficients<Lexicographical,T,Index2D> &coeff, const Basis &basis2d, const char* filename);

template <typename T>
void
readCoeffVector2D(Coefficients<Lexicographical,T,Index2D>& coeff, const char* filename, bool append = false);

template <typename T, typename Basis>
void
saveIndexSet2D(const IndexSet<Index2D> &indexset, const Basis &basis2d, const char* filename);
  
template <typename T>
void
readIndexSet2D(IndexSet<Index2D>& indexset, const char* filename, bool append = false);

} // namespace lawa

#include <lawa/methods/rb/postprocessing/plotting.tcc>

#endif // LAWA_METHODS_RB_POSTPROCESSING_PLOTTING_H
