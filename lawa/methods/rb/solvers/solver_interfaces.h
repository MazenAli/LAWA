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

#ifndef LAWA_METHODS_RB_SOLVERS_SOLVER_INTERFACES_H_
#define LAWA_METHODS_RB_SOLVERS_SOLVER_INTERFACES_H_

namespace lawa {

template<typename DataType, typename LHS, typename RHS,
typename TrialBasis, typename TestBasis, typename Params>
struct GeneralSolverInterface {

    typedef LHS 						LHSType;
    typedef RHS							RHSType;
    typedef TrialBasis					TrialBasisType;
    typedef TestBasis					TestBasisType;
    typedef Params						ParamType;

	/* Returns a solution vector u
	 * 	(as usual still with preconditioner)
	 */
	DataType
	solve() = 0;

	/* Undo the preconditioning of (solution) vector u
	 *  (in PG settings, we assume u \in TrialSpace)
	 */
	void
	remove_preconditioner(DataType& u) = 0;

	/* Get (affine) left/right hand side
	 * in order to be able to set parameters
	 */
	LHS&
	get_lhs() = 0;

	RHS&
	get_rhs() = 0;

    const TrialBasis&
    get_trialbasis();

    const TestBasis&
    get_testbasis();

    Params	params;
};

} // namespace lawa

#endif /* LAWA_METHODS_RB_SOLVERS_SOLVER_INTERFACES_H_ */
