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
 
#ifndef LAWA_METHODS_RB_POSTPROCESSING_ERROR_FCTS_H
#define LAWA_METHODS_RB_POSTPROCESSING_ERROR_FCTS_H

namespace lawa {

template<typename T, typename Basis2D, typename Prec>
T
L2_H1_error(Basis2D& basis2d, const Coefficients<Lexicographical,T,Index2D> &u, Prec& P,
            T (*u_ref)(T,T), T (*dx_u_ref)(T,T), T a_t, T b_t, int n_t, T a_x, T b_x, int n_x);

} // namespace lawa

#include <lawa/methods/rb/postprocessing/error_fcts.tcc>

#endif // LAWA_METHODS_RB_POSTPROCESSING_ERROR_FCTS_H
