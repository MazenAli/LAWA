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
 
#ifndef LAWA_METHODS_UNIFORM_POSTPROCESSING_EVALUATE_H
#define LAWA_METHODS_UNIFORM_POSTPROCESSING_EVALUATE_H 1

namespace lawa {

template <typename X, typename Basis>
typename X::ElementType
evaluate(const Basis &basis, const int J_x, const int J_y, const flens::DenseVector<X> &coeffs,
         const typename X::ElementType x, const typename X::ElementType y, const int deriv_x,
         const int deriv_y);

} // namespace lawa

#include <lawa/methods/uniform/postprocessing/evaluate.tcc>


#endif // LAWA_METHODS_UNIFORM_POSTPROCESSING_EVALUATE_H

