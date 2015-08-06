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

#ifndef  LAWA_METHODS_ADAPTIVE_SOLVERS_MULTITREEAWGM_PG_H
#define  LAWA_METHODS_ADAPTIVE_SOLVERS_MULTITREEAWGM_PG_H 1

#include <lawa/methods/adaptive/datastructures/datastructures.h>
#include <lawa/methods/adaptive/solvers/solver_parameters.h>

namespace lawa {

/**
 * Class for the solution of problems using an Adaptive Wavelet Galerkin Method
 * based on multitree matrix-vector-multiplications
 * Can handle different trial and test bases as well as left/right preconditioning
 *
 * !!! Assumes that neither the Operators nor the right hand side are already
 * !!! preconditioned.
 */
template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
class MultiTreeAWGM_PG {

    typedef typename LocalOperator::T 	T;
    typedef T (*sol_fct_2d)(T,T);

public:

    typedef LocalOperator 				LHSType;
    typedef RHS							RHSType;
    typedef TrialBasis					TrialBasisType;
    typedef TestBasis					TestBasisType;
    typedef AWGM_PG_Parameters			ParamType;
    typedef TrialPrec					TrialPrecType;
    typedef TestPrec					TestPrecType;

    MultiTreeAWGM_PG(const TrialBasis &_trialbasis, const TestBasis& _testbasis,
    				LocalOperator &_Op, LocalOperatorTransp& _OpTransp, RHS &_F,
    				TrialPrec &_trialPrec, TestPrec& _testPrec);

    MultiTreeAWGM_PG(const TrialBasis &_trialbasis, const TestBasis& _testbasis,
    				LocalOperator &_Op, LocalOperatorTransp& _OpTransp, RHS &_F,
    				TrialPrec &_trialPrec, TestPrec& _testPrec,
    				AWGM_PG_Parameters& _awgm_params, IS_Parameters& _is_params);

    // CGLS solve
    void
    solve(Coefficients<Lexicographical,T,Index> &u, IndexSet<Index>& init_Lambda_trial, IndexSet<Index>& init_Lambda_test);

    // CGLS solve
    void
    solve(Coefficients<Lexicographical,T,Index> &u);

    // CGLS solve
    Coefficients<Lexicographical,T,Index>
    solve();

    RHS&
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
    set_initial_indexsets(const IndexSet<Index> _LambdaTrial, const IndexSet<Index> _LambdaTest);

    void
    remove_preconditioner(Coefficients<Lexicographical,T,Index> &u);

    ParamType&
    access_params();

    void
    reset_info();

    AWGM_PG_Parameters						awgm_params;
    IS_Parameters							is_params;

private:

    const TrialBasis&                       trialbasis;
    const TestBasis&                        testbasis;
    LocalOperator&                          Op;
    LocalOperatorTransp&					OpTransp;
    RHS&                                    F;
    TrialPrec&                          	trialPrec;
    TestPrec&                          		testPrec;

    AWGM_PG_Information						awgm_info;

    sol_fct_2d								exact_sol;

    IndexSet<Index> 						default_init_Lambda_trial,
    										default_init_Lambda_test;

    MultiTreeAWGM_PG(const MultiTreeAWGM_PG&);
};


}    //namespace lawa

#include <lawa/methods/adaptive/solvers/multitreeawgm_pg.tcc>

#endif    // LAWA_METHODS_ADAPTIVE_SOLVERS_MULTITREEAWGM_PG_H

