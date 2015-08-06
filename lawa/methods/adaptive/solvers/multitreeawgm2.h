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

#ifndef  LAWA_METHODS_ADAPTIVE_SOLVERS_MULTITREEAWGM2_H
#define  LAWA_METHODS_ADAPTIVE_SOLVERS_MULTITREEAWGM2_H 1

#include <lawa/methods/adaptive/datastructures/datastructures.h>
#include <lawa/methods/adaptive/solvers/solver_parameters.h>

namespace lawa {

/**
 * Class for the solution of problems using an Adaptive Wavelet Galerkin Method
 * based on multitree matrix-vector-multiplications
 *
 * !!! Assumes that neither the Operators nor the right hand side are already
 * !!! preconditioned.
 */
template <typename Index, typename Basis, typename LocalOperator,
		  typename RHS, typename Preconditioner>
class MultiTreeAWGM2 {

    typedef typename LocalOperator::T 	T;
    typedef AWGM_Parameters				ParamType;

    typedef T (*sol_fct_2d)(T,T);

public:

    typedef LocalOperator 			  	LHSType;
    typedef RHS						 	RHSType;
    typedef Basis						BasisType;
    typedef Basis						TrialBasisType;
    typedef Basis						TestBasisType;
    typedef Preconditioner				PrecType;
    typedef Preconditioner				TrialPrecType;
    typedef Preconditioner				TestPrecType;

    MultiTreeAWGM2(const Basis &_basis, LocalOperator &_Op, RHS &_F, Preconditioner &_Prec);

    MultiTreeAWGM2(const Basis &_basis, LocalOperator &_Op, RHS &_F, Preconditioner &_Prec,
    				AWGM_Parameters& _awgm_params, IS_Parameters& _is_params);

    // CG solve (returns res_norm)
    T
    solve(Coefficients<Lexicographical,T,Index> &u, IndexSet<Index>& init_Lambda, T old_res_norm = 1.);

    // CG solve
    T
    solve(Coefficients<Lexicographical,T,Index> &u, T old_res_norm = 1.);

    // CG solve
    Coefficients<Lexicographical,T,Index>
    solve();

    void
    coarsen(Coefficients<Lexicographical,T,Index>& u, T delta);

    void
    coarsen(Coefficients<Lexicographical,T,Index>& u);

    void
    set_sol(sol_fct_2d _sol);

    void
    set_initial_indexset(const IndexSet<Index> _LambdaTrial);

    void
    remove_preconditioner(Coefficients<Lexicographical,T,Index> &u);

    T
    estimate_residual_norm(Coefficients<Lexicographical,T,Index>& u);

    RHS&
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

    AWGM_Parameters						awgm_params;
    IS_Parameters						is_params;

private:

    const Basis&                       basis;
    LocalOperator&                     Op;
    RHS&                               F;
    Preconditioner&                    Prec;

    AWGM_Information				   awgm_info;

    sol_fct_2d						   exact_sol;

    IndexSet<Index>					   default_init_Lambda;

    MultiTreeAWGM2(const MultiTreeAWGM2&);
};


}    //namespace lawa

#include <lawa/methods/adaptive/solvers/multitreeawgm2.tcc>

#endif    // LAWA_METHODS_ADAPTIVE_SOLVERS_MULTITREEAWGM2_H

