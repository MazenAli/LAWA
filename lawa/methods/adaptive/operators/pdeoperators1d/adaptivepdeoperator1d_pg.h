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


#ifndef LAWA_METHODS_ADAPTIVE_OPERATORS_ADAPTIVEPDEOPERATOR1D_PG_H
#define LAWA_METHODS_ADAPTIVE_OPERATORS_ADAPTIVEPDEOPERATOR1D_PG_H 1

#include <lawa/settings/enum.h>
#include <lawa/settings/typetraits.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/operators/pdeoperators1d/pdeoperator1d_pg.h>
#include <lawa/preconditioners/preconditioners.h>

namespace lawa {

/* Adaptive PDE OPERATOR: Petrov Galerkin version
 *
 *    a(v,u) =  diffusion * Integral(v_x * u_x) +  convection * Integral(v * u_x)
 *              + reaction * Integral(v * u)
 *
 */
template <typename T, typename TrialBasis, typename TestBasis>
struct AdaptivePDEOperator1D_PG
{
    typedef flens::SparseGeMatrix<flens::extensions::CRS<T,flens::CRS_General> >        SparseMatrixT;
    typedef NoCompression<T,Index1D,TrialBasis>                        		NoCompression1D;
    typedef NoPreconditioner<T,Index1D>                                     NoPreconditioner1D;
    typedef PDEOperator1D_PG<T,TrialBasis,TestBasis>                   		PDEOp1D;
    typedef MapMatrix<T, Index1D, PDEOp1D,
                     NoCompression1D, NoPreconditioner1D>                   DataIdentityOp1D;


    AdaptivePDEOperator1D_PG(const TrialBasis& _trialbasis, const TestBasis& _testbasis,
    							  T _reaction, T _convection, T _diffusion);

    T
    operator()(const Index1D &row_index, const Index1D &col_index);

    T
    operator()(XType xtype_row, FLENS_DEFAULT_INDEXTYPE j_row, FLENS_DEFAULT_INDEXTYPE k_row, XType xtype_col, FLENS_DEFAULT_INDEXTYPE j_col, FLENS_DEFAULT_INDEXTYPE k_col);

    void
    toFlensSparseMatrix(const IndexSet<Index1D>& LambdaRow, const IndexSet<Index1D>& LambdaCol,
                        SparseMatrixT &A_flens, FLENS_DEFAULT_INDEXTYPE J=-1, bool useLinearIndex=false);

    void
    clear();

    const TrialBasis&           trialbasis1d;
    const TestBasis&            testbasis1d;
    NoCompression1D             compression1d;
    const PDEOp1D          		pde_op1d;
    NoPreconditioner1D          prec1d;
    DataIdentityOp1D            data;

};

}   //namespace lawa

#include <lawa/methods/adaptive/operators/pdeoperators1d/adaptivepdeoperator1d_pg.tcc>

#endif // LAWA_METHODS_ADAPTIVE_OPERATORS_ADAPTIVEPDEOPERATOR1D_PG_H

