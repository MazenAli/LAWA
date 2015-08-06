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

#ifndef LAWA_METHODS_ADAPTIVE_ALGORITHMS_INDEXSET_GENERATION_H_
#define LAWA_METHODS_ADAPTIVE_ALGORITHMS_INDEXSET_GENERATION_H_

#include <lawa/methods/adaptive/datastructures/datastructures.h>

namespace lawa {

template <typename T, typename Basis2D>
void
getSparseGridIndexSet(const Basis2D &basis, IndexSet<Index2D> &Lambda, int j, int deltaL, T gamma = 0.);


template <typename Basis1D>
void
getFullIndexSet(const Basis1D &basis, IndexSet<Index1D> &Lambda, int J);


template <typename Basis2D>
void
getFullIndexSet(const Basis2D &basis, IndexSet<Index2D> &Lambda, int J1, int J2, int deltaL);

template <typename Basis2D>
void
getScalingFctIndexSet(const Basis2D &basis, IndexSet<Index2D> &Lambda, int J1, int J2);

} // namespace lawa

#include <lawa/methods/adaptive/algorithms/indexset_generation.tcc>

#endif /* LAWA_METHODS_ADAPTIVE_ALGORITHMS_INDEXSET_GENERATION_H_ */
