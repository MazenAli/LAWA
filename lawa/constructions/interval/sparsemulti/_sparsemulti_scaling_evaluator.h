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
 
#ifndef LAWA_CONSTRUCTIONS_INTERVAL_SPARSEMULTI_SCALING_EVALUATOR_H
#define LAWA_CONSTRUCTIONS_INTERVAL_SPARSEMULTI_SCALING_EVALUATOR_H 1

namespace lawa {

//--- cubic evaluators -----------------------------------------------------

template <typename T>
    T
    _cubic_sparsemulti_scaling_left_evaluator0(T x, unsigned short deriv);

template <typename T>
    T
    _cubic_sparsemulti_scaling_inner_evaluator0(T x, unsigned short deriv);

template <typename T>
    T
    _cubic_sparsemulti_scaling_inner_evaluator1(T x, unsigned short deriv);

template <typename T>
    T
    _cubic_sparsemulti_scaling_right_evaluator0(T x, unsigned short deriv);

} // namespace lawa

#include <lawa/constructions/interval/sparsemulti/_sparsemulti_scaling_evaluator.tcc>

#endif // LAWA_CONSTRUCTIONS_INTERVAL_SPARSEMULTI_SCALING_EVALUATOR_H
