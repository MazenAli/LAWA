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

#ifndef LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_FLEXIBLECOMPOUNDLOCALOPERATOR_H
#define LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_FLEXIBLECOMPOUNDLOCALOPERATOR_H 1

#include <vector>
#include <lawa/methods/adaptive/datastructures/datastructures.h>

namespace lawa {

template <typename Index, typename LocalOperatorType>
class FlexibleCompoundLocalOperator {

public:

	typedef typename LocalOperatorType::T T;

	FlexibleCompoundLocalOperator(std::vector<LocalOperatorType*>& _localops);

	void
	eval(const Coefficients<Lexicographical,T,Index> &v,
		 Coefficients<Lexicographical,T,Index> &Av);

	template <typename Preconditioner>
	void
	eval(Coefficients<Lexicographical,T,Index> &v,
		 Coefficients<Lexicographical,T,Index> &Av, Preconditioner &P);

	template <typename RightPrec, typename LeftPrec>
	void
	eval(Coefficients<Lexicographical,T,Index> &v,
		 Coefficients<Lexicographical,T,Index> &Av, RightPrec& rightP, LeftPrec& leftP);

	void
	eval(Coefficients<Lexicographical,T,Index> &v,
		 Coefficients<Lexicographical,T,Index> &Av,
		 Coefficients<Lexicographical,T,Index>& rightP,
		 Coefficients<Lexicographical,T,Index>& leftP);

	void
	eval(std::size_t i,
		 const Coefficients<Lexicographical,T,Index> &v,
		 Coefficients<Lexicographical,T,Index> &Av);

	T
	eval(std::size_t i, const Index& ind_row, const Index& ind_col);

	std::size_t
	size();
	
	void
    clear();

protected:

	std::vector<LocalOperatorType*>& localops;

private:

	typedef typename Coefficients<Lexicographical,T,Index>::iterator    coeff_it;
    typedef typename Coefficients<Lexicographical,T,Index>::iterator    const_coeff_it;
    typedef typename IndexSet<Index1D>::const_iterator                  const_set1d_it;


	FlexibleCompoundLocalOperator(const FlexibleCompoundLocalOperator& rhs);

};

}   // namespace lawa

#include <lawa/methods/adaptive/operators/localoperators/flexiblecompoundlocaloperator.tcc>

#endif // LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_FLEXIBLECOMPOUNDLOCALOPERATOR_H
