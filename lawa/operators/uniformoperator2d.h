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
 
#ifndef LAWA_OPERATORS_UNIFORMOPERATOR2D_H
#define LAWA_OPERATORS_UNIFORMOPERATOR2D_H 1

#include <lawa/operators/operator2d.h>

namespace lawa {
    
template <typename T>
struct UniformOperator2D : public Operator2D<T> {
	
    virtual T
    operator()(XType row_xtype_x, int j1_x, long k1_x,
               XType row_xtype_y, int j1_y, long k1_y,
               XType col_xtype_x, int j2_x, long k2_x,
               XType col_xtpye_y, int j2_y, long k2_y) const = 0;
    
    virtual T
    operator()(const Index2D &row_index, const Index2D &col_index) const = 0;
    
};
    
} // namespace lawa

#endif // LAWA_OPERATORS_UNIFORMOPERATOR2D_H
