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

#ifndef LAWA_METHODS_RB_RIGHTHANDSIDES_AFFINERHS_H_
#define LAWA_METHODS_RB_RIGHTHANDSIDES_AFFINERHS_H_

#include <lawa/methods/adaptive/righthandsides/flexiblecompoundrhs.h>

namespace lawa {

template <typename T, typename Index, typename RHSType, typename ParamType>
class AffineRhs : public FlexibleCompoundRhs<T, Index,RHSType>{

public:
	AffineRhs(ThetaStructure<ParamType>& _thetas,
			  std::vector<RHSType*>& _rhsvec);

	T
	operator()(const Index& index);

	Coefficients<Lexicographical,T,Index>
	operator()(const IndexSet<Index> &indexset);

	void
	set_param(ParamType& mu);

private:

	ThetaStructure<ParamType>& thetas;

	AffineRhs(const AffineRhs& rhs);
};

} // namespace lawa

#include <lawa/methods/rb/righthandsides/affinerhs.tcc>

#endif /* LAWA_METHODS_RB_RIGHTHANDSIDES_AFFINERHS_H_ */
