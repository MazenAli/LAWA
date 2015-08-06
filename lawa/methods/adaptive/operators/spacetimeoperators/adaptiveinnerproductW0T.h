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

#ifndef LAWA_METHODS_ADAPTIVE_OPERATORS_ADAPTIVEINNERPRODUCTW0T_PG_H
#define LAWA_METHODS_ADAPTIVE_OPERATORS_ADAPTIVEINNERPRODUCTW0T_PG_H 1

#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/compressions/nocompression.h>
#include <lawa/methods/adaptive/compressions/compression_pde1d.h>
#include <lawa/methods/adaptive/compressions/compression_pde2d.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/mapmatrix.h>
#include <lawa/methods/adaptive/operators/pdeoperators2d/adaptiveoperator2d.h>
#include <lawa/operators/pdeoperators1d/identityoperator1d_pg.h>
#include <lawa/operators/pdeoperators1d/laplaceoperator1d_pg.h>
#include <lawa/operators/pdeoperators1d/convectionoperator1d_pg.h>
#include <lawa/operators/spacetimeoperators/spacetimeoperators.h>
#include <lawa/preconditioners/nopreconditioner.h>
#include <lawa/preconditioners/spacetimepreconditioners/spacetimepreconditioners.h>

namespace lawa {

/* Space-Time Innner Product in W(0,T) : Petrov Galerkin version
 *
 *  a(v,u) =  Integral(v1 * u1) * [Integral(v2 * u2) + Integral(v2_x * u2_x)]
 *          + Integral(v1_t * u1_t) * [Integral(v2_x * u2_x)]_V'
 *
 *  Template Parameters:
 *      TrialPrec :        	left preconditioner
 *      TestPrec:        	right preconditioner
 */
template <typename T, typename Basis2D, typename TrialPrec, typename TestPrec>
struct AdaptiveInnerProductW0T : public AdaptiveOperator2D<T> {

	typedef flens::SparseGeMatrix<flens::extensions::CRS<T,flens::CRS_General> >                  SparseMatrixT;

    typedef typename Basis2D::FirstBasisType    Basis_t;
    typedef typename Basis2D::SecondBasisType   Basis_x;

    /*typedef CompressionPDE1D<T, TrialBasis_t>      Compression1D_Trial_t;
    typedef CompressionPDE1D<T, TrialBasis_x>      Compression1D_Trial_x;
    typedef CompressionPDE1D<T, TestBasis_t>       Compression1D_Test_t;
    typedef CompressionPDE1D<T, TestBasis_x>       Compression1D_Test_x;
    typedef CompressionPDE2D<T, TrialBasis>        Compression2D_Trial;
    typedef CompressionPDE2D<T, TestBasis>         Compression2D_Test;
    */
    typedef NoCompression<T, Index1D, Basis_t>   Compression1D_t;
    typedef NoCompression<T, Index1D, Basis_x>   Compression1D_x;
    typedef NoCompression<T, Index2D, Basis2D>     Compression2D;

    typedef NoPreconditioner<T,Index1D>         NoPreconditioner1D;

    typedef IdentityOperator1D<T, Basis_t>      IdentityOperator_t;
    typedef IdentityOperator1D<T, Basis_x>      IdentityOperator_x;
    typedef LaplaceOperator1D<T, Basis_t>       LaplaceOperator_t;
    typedef LaplaceOperator1D<T, Basis_x>       LaplaceOperator_x;

    typedef MapMatrix<T, Index1D, IdentityOperator_t,
                               Compression1D_t, NoPreconditioner1D>   DataIdentity_t;
    typedef MapMatrix<T, Index1D, IdentityOperator_x,
                               Compression1D_x, NoPreconditioner1D>   DataIdentity_x;
    typedef MapMatrix<T, Index1D, LaplaceOperator_t,
                               Compression1D_t, NoPreconditioner1D>   DataLaplace_t;
    typedef MapMatrix<T, Index1D, LaplaceOperator_x,
                               Compression1D_x, NoPreconditioner1D>   DataLaplace_x;

    AdaptiveInnerProductW0T(const Basis2D& _basis, TrialPrec& _trialprec, TestPrec& _testprec);

    // call of p_left * a_operator * p_right
    T
    operator()(const Index2D &row_index, const Index2D &col_index);

    virtual Coefficients<Lexicographical,T,Index2D>
    mv(const IndexSet<Index2D> &LambdaRow,
       const Coefficients<Lexicographical,T,Index2D> &x);

	void
	toFlensSparseMatrix(const IndexSet<Index2D> &LambdaRow, const IndexSet<Index2D> &LambdaCol,
						SparseMatrixT &A, T tol, bool useLinearIndex=false);

    void
    clear();

    const Basis2D&   basis;

    Compression1D_t     compression_1d_t;
    Compression1D_x     compression_1d_x;
    Compression2D       compression;

    Coefficients<Lexicographical,T,Index2D> P_trial_data;
    Coefficients<Lexicographical,T,Index2D> P_test_data;

    const TrialPrec&   trialprec;
    const TestPrec&    testprec;
    NoPreconditioner1D  noprec;

    const IdentityOperator_t    op_identity_t;
    const IdentityOperator_x    op_identity_x;
    const LaplaceOperator_t     op_laplace_t;
    const LaplaceOperator_x     op_laplace_x;

    DataIdentity_t      data_identity_t;
    DataIdentity_x      data_identity_x;
    DataLaplace_t       data_laplace_t;
    DataLaplace_x       data_laplace_x;

};

} // namespace lawa

#include <lawa/methods/adaptive/operators/spacetimeoperators/adaptiveinnerproductW0T.tcc>


#endif /* LAWA_METHODS_ADAPTIVE_OPERATORS_ADAPTIVEINNERPRODUCTW0T_PG_H */
