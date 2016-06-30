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

#ifndef  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_MATRIXOPERATIONS_H
#define  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_MATRIXOPERATIONS_H 1

#include <extensions/flens/cg.h>
#include <extensions/flens/gmres.h>
#include <lawa/flensforlawa.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>

namespace lawa {


template <typename T, typename Index, typename MA>
void
toFlensSparseMatrix(MA &A, const IndexSet<Index>& LambdaRow, const IndexSet<Index>& LambdaCol,
                    flens::SparseGeMatrix<flens::extensions::CRS<T,flens::CRS_General> > &A_flens);

template <typename T, typename Index, typename SpaceIndex, typename MA>
void
toFlensSparseMatrix(MA &A, const IndexSet<Index>& LambdaRowOp, const IndexSet<SpaceIndex>& LambdaRowInitCond,
                    const IndexSet<Index>& LambdaCol,
                    flens::SparseGeMatrix<flens::extensions::CRS<T,flens::CRS_General> > &A_flens);

template <typename T, typename RowIndex, typename ColIndex, typename MA>
Coefficients<Lexicographical,T,RowIndex>
mv(const IndexSet<RowIndex> &LambdaRow, MA &A, const Coefficients<Lexicographical,T,ColIndex > &v);

//requires lambdaTilde!!!
template <typename T, typename Index, typename MA>
Coefficients<Lexicographical,T,Index>
mv_sparse(const IndexSet<Index> &LambdaRow, MA &A, const Coefficients<Lexicographical,T,Index > &v);

/*  Attention, time t is not checked for storing entries!!!
template <typename T, typename Index, typename MA>
Coefficients<Lexicographical,T,Index>
mv(T t, const IndexSet<Index> &LambdaRow, MA &A, const Coefficients<Lexicographical,T,Index > &v);

template <typename T, typename MA>
Coefficients<Lexicographical,T,Index2D>
mv_sparse(t, const IndexSet<Index2D> &LambdaRow, MA &A, const Coefficients<Lexicographical,T,Index2D > &v);
*/

/*
template <typename T, typename Index, typename MA>
FLENS_DEFAULT_INDEXTYPE
CG_Solve(const IndexSet<Index> &Lambda, MA &A, Coefficients<Lexicographical,T,Index > &u,
         const Coefficients<Lexicographical,T,Index > &f, T &res, T tol = 1e-6,
         bool useOptimizedAssembling=false, FLENS_DEFAULT_INDEXTYPE maxIterations = 1000);

template <typename T, typename Index, typename MA>
FLENS_DEFAULT_INDEXTYPE
GMRES_Solve(const IndexSet<Index> &Lambda, MA &A, Coefficients<Lexicographical,T,Index > &u, 
            const Coefficients<Lexicographical,T,Index > &f, T &res, T tol = 10e-6,
            bool useOptimizedAssembling=false, FLENS_DEFAULT_INDEXTYPE maxIterations = 1000);
            
template <typename T, typename Index, typename MA>
FLENS_DEFAULT_INDEXTYPE
GMRESM_Solve(const IndexSet<Index> &Lambda, MA &A, Coefficients<Lexicographical,T,Index > &u,
             const Coefficients<Lexicographical,T,Index > &f, T &residual, T tol=10e-6,
             bool useOptimizedAssembling=false, FLENS_DEFAULT_INDEXTYPE maxIterations=1000, FLENS_DEFAULT_INDEXTYPE m=20);

template <typename T, typename Index, typename MA>
FLENS_DEFAULT_INDEXTYPE
GMRES_Solve_PG(const IndexSet<Index> &Lambda, MA &A, Coefficients<Lexicographical,T,Index > &u, 
            const Coefficients<Lexicographical,T,Index > &f, T &res, T tol = 10e-6, FLENS_DEFAULT_INDEXTYPE maxIterations = 1000);

template <typename T, typename Index, typename SpaceIndex, typename MA>
FLENS_DEFAULT_INDEXTYPE
CGLS_Solve(const IndexSet<Index> &LambdaRowOp,
           const IndexSet<SpaceIndex> &LambdaRowInitCond, MA &A,
           const IndexSet<Index> &LambdaCol,
           Coefficients<Lexicographical,T,Index > &u,
           const Coefficients<Lexicographical,T,Index > &f,
           const Coefficients<Lexicographical,T,SpaceIndex > &u0,
           T &res, T tol = 1e-6, FLENS_DEFAULT_INDEXTYPE maxIterations = 1000);

template <typename T, typename Index, typename MA>
FLENS_DEFAULT_INDEXTYPE
CGLS_Solve(const IndexSet<Index> &LambdaRow, const IndexSet<Index> &LambdaCol,
           MA &A, Coefficients<Lexicographical,T,Index > &u,
           const Coefficients<Lexicographical,T,Index > &f,
           T &res, T tol = 1e-6, FLENS_DEFAULT_INDEXTYPE maxIterations = 1000);
*/
} // namespace lawa

#include <lawa/methods/adaptive/datastructures/matrixoperations.tcc>

#endif    // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_MATRIXOPERATIONS_H

