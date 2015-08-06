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

#ifndef LAWA_METHODS_ADAPTIVE_OPERATORS_ADAPTIVEOPERATOR2D_H
#define LAWA_METHODS_ADAPTIVE_OPERATORS_ADAPTIVEOPERATOR2D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/operators/operator2d.h>

namespace lawa {
  
  template <typename T>
  struct AdaptiveOperator2D : Operator2D<T> {
    
    typedef flens::SparseGeMatrix<flens::extensions::CRS<T,flens::CRS_General> >                  SparseMatrixT;

    virtual T
    operator()(const Index2D &row_index, const Index2D &col_index) = 0;
    
    virtual Coefficients<Lexicographical,T,Index2D>
    mv(const IndexSet<Index2D> &LambdaRow,
       const Coefficients<Lexicographical,T,Index2D> &x) = 0;
    
    virtual void
    toFlensSparseMatrix(const IndexSet<Index2D> &LambdaRow, const IndexSet<Index2D> &LambdaCol, 
												SparseMatrixT &A, T eps, bool useLinearIndex=false) = 0; 
    
  };
  
} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_OPERATORS_ADAPTIVEOPERATOR2D_H
