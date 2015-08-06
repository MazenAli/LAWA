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


#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVELAPLACEOPERATOR1D_PG_H
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVELAPLACEOPERATOR1D_PG_H 1

#include <lawa/settings/enum.h>
#include <lawa/settings/typetraits.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/mapmatrix.h>
#include <lawa/operators/pdeoperators1d/laplaceoperator1d_pg.h>
#include <lawa/preconditioners/preconditioners.h>


namespace lawa {

/* Adaptive Laplace OPERATOR: Petrov Galerkin version
 *
 *    a(v,u) = Integral(v_x * u_x)
 *
 */
template <typename T, typename TrialBasis, typename TestBasis>
struct AdaptiveLaplaceOperator1D_PG
{
    typedef flens::SparseGeMatrix<flens::extensions::CRS<T,flens::CRS_General> >        SparseMatrixT;
    typedef NoCompression<T,Index1D,TrialBasis>                        		NoCompression1D;
    typedef NoPreconditioner<T,Index1D>                                     NoPreconditioner1D;
    typedef LaplaceOperator1D_PG<T, TrialBasis,TestBasis>                   LaplaceOp1D;
    typedef MapMatrix<T, Index1D, LaplaceOp1D,
                     NoCompression1D, NoPreconditioner1D>              DataLaplaceOp1D;


    AdaptiveLaplaceOperator1D_PG(const TrialBasis& _trialbasis, const TestBasis& _testbasis);

    T
    operator()(const Index1D &row_index, const Index1D &col_index);

    T
    operator()(XType xtype_row, int j_row, long k_row, XType xtype_col, int j_col, long k_col);

    void
    clear();

    const TrialBasis&           trialbasis1d;
    const TestBasis&            testbasis1d;
    NoCompression1D             compression1d;
    const LaplaceOp1D          laplace_op1d;
    NoPreconditioner1D         prec1d;
    DataLaplaceOp1D            data;

};

}   //namespace lawa

#include <lawa/methods/adaptive/operators/pdeoperators1d/adaptivelaplaceoperator1d_pg.tcc>

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVELAPLACEOPERATOR1D_PG_H

