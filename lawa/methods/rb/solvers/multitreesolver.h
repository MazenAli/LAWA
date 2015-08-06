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

#ifndef  LAWA_METHODS_RB_SOLVERS_MULTITREESOLVER_H
#define  LAWA_METHODS_RB_SOLVERS_MULTITREESOLVER_H 1

#include <lawa/methods/adaptive/datastructures/datastructures.h>
#include <lawa/methods/adaptive/solvers/solver_parameters.h>

namespace lawa {

/**
* Class for the solution of problems using an Wavelet Galerkin Method
* on FIXED index sets,
* based on multitree matrix-vector-multiplications
 *
 * !!! Assumes that neither the Operators nor the right hand side are already
 * !!! preconditioned.
 */
template <typename Index, typename Basis, typename LocalOperator,
		  typename F_RHS, typename Preconditioner>
class MultiTreeSolver {

    typedef typename LocalOperator::T 	T;
    typedef LocalOperator 			  	LHSType;
    typedef F_RHS						RHSType;
    typedef Basis						BasisType;
    typedef Basis						TrialBasisType;
    typedef Basis						TestBasisType;
    typedef Preconditioner				PrecType;
    typedef Preconditioner				TrialPrecType;
    typedef Preconditioner				TestPrecType;
    typedef ISWGM_Parameters			ParamType;



    typedef T (*sol_fct_2d)(T,T);

public:

    MultiTreeSolver(const Basis &_basis, LocalOperator &_Op, F_RHS &_F, Preconditioner &_Prec);

    MultiTreeSolver(const Basis &_basis, LocalOperator &_Op, F_RHS &_F, Preconditioner &_Prec,
    				ISWGM_Parameters& _iswgm_params, IS_Parameters& _is_params);

    // CG solve
    T
    solve(Coefficients<Lexicographical,T,Index> &u, IndexSet<Index>& init_Lambda, T dummy = 0.);

    // CG solve
    T
    solve(Coefficients<Lexicographical,T,Index> &u, T dummy = 0.);

    // CG solve
    Coefficients<Lexicographical,T,Index>
    solve();

    void
    set_sol(sol_fct_2d _sol);

    void
    set_indexset(const IndexSet<Index> _LambdaTrial);

    void
    remove_preconditioner(Coefficients<Lexicographical,T,Index> &u);

    F_RHS&
    get_rhs();

    LocalOperator&
    get_lhs();

    const Basis&
    get_trialbasis();

    const Basis&
    get_testbasis();

    Preconditioner&
    get_trialprec();

    Preconditioner&
    get_testprec();

    ParamType&
    access_params();

    void
    reset_info();

    ISWGM_Parameters				   iswgm_params;

    IS_Parameters					   is_params;
    
    void
    coarsen(Coefficients<Lexicographical,T,Index>& u);

private:

    const Basis&                       basis;
    LocalOperator&                     Op;
    F_RHS&                             F;
    Preconditioner&                    Prec;

    sol_fct_2d						   exact_sol;

    ISWGM_Information				   iswgm_info;

    IndexSet<Index>					   default_Lambda;

    MultiTreeSolver(const MultiTreeSolver&);
};


}    //namespace lawa

#include <lawa/methods/rb/solvers/multitreesolver.tcc>

#endif    // LAWA_METHODS_RB_SOLVERS_MULTITREESOLVER_H

