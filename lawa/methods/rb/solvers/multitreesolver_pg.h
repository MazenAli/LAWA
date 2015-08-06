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

#ifndef  LAWA_METHODS_RB_SOLVERS_MULTITREESOLVER_PG_H
#define  LAWA_METHODS_RB_SOLVERS_MULTITREESOLVER_PG_H 1

#include <lawa/methods/adaptive/datastructures/datastructures.h>
#include <lawa/methods/adaptive/solvers/solver_parameters.h>

namespace lawa {

/**
 * Class for the solution of problems using an Wavelet Galerkin Method
 * on FIXED index sets,
 * based on multitree matrix-vector-multiplications
 * Can handle different trial and test bases as well as left/right preconditioning
 *
 * !!! Assumes that neither the Operators nor the right hand side are already
 * !!! preconditioned.
 */
template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename F_RHS,
		  typename TrialPrec, typename TestPrec>
class MultiTreeSolver_PG {

    typedef typename LocalOperator::T 	T;
    typedef T (*sol_fct_2d)(T,T);

public:

    typedef LocalOperator 				LHSType;
    typedef F_RHS						RHSType;
    typedef TrialBasis					TrialBasisType;
    typedef TestBasis					TestBasisType;
    typedef TrialPrec   				TrialPrecType;
    typedef TestPrec    				TestPrecType;
    typedef ISWGM_Parameters			ParamType;

    MultiTreeSolver_PG(const TrialBasis &_trialbasis, const TestBasis& _testbasis,
    				LocalOperator &_Op, LocalOperatorTransp& _OpTransp, F_RHS &_F,
    				TrialPrec &_trialPrec, TestPrec& _testPrec);

    MultiTreeSolver_PG(const TrialBasis &_trialbasis, const TestBasis& _testbasis,
    				LocalOperator &_Op, LocalOperatorTransp& _OpTransp, F_RHS &_F,
    				TrialPrec &_trialPrec, TestPrec& _testPrec,
    				ISWGM_Parameters& _iswgm_params, IS_Parameters& _is_params);

    // CGLS solve
    T
    solve(Coefficients<Lexicographical,T,Index> &u, IndexSet<Index>& Lambda_trial, IndexSet<Index>& Lambda_test, T dummy = 0.);

    // CGLS solve
    T
    solve(Coefficients<Lexicographical,T,Index> &u, T dummy = 0.);

    // CGLS solve
    Coefficients<Lexicographical,T,Index>
    solve();

    F_RHS&
    get_rhs();

    LocalOperator&
    get_lhs();

    const TrialBasis&
    get_trialbasis();

    const TestBasis&
    get_testbasis();
    
    TrialPrec&
    get_trialprec();

    TestPrec&
    get_testprec();

    void
    set_sol(sol_fct_2d _sol);

    void
    set_indexsets(const IndexSet<Index> _LambdaTrial, const IndexSet<Index> _LambdaTest);

    void
    remove_preconditioner(Coefficients<Lexicographical,T,Index> &u);

    ParamType&
    access_params();

    void
    reset_info();

    ISWGM_Parameters						iswgm_params;
    IS_Parameters							is_params;
    
    void
    coarsen(Coefficients<Lexicographical,T,Index>& u);

private:

    const TrialBasis&                       trialbasis;
    const TestBasis&                        testbasis;
    LocalOperator&                          Op;
    LocalOperatorTransp&					OpTransp;
    F_RHS&                                  F;
    TrialPrec&                          	trialPrec;
    TestPrec&                          		testPrec;

    sol_fct_2d								exact_sol;

    ISWGM_PG_Information					iswgm_info;

    IndexSet<Index> 						default_Lambda_trial,
    										default_Lambda_test;

    MultiTreeSolver_PG(const MultiTreeSolver_PG&);
};


}    //namespace lawa

#include <lawa/methods/rb/solvers/multitreesolver_pg.tcc>

#endif    // LAWA_METHODS_RB_SOLVERS_MULTITREESOLVER_PG_H

