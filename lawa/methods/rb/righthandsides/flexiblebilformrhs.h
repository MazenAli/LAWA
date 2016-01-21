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

#ifndef LAWA_METHODS_RB_RIGHTHANDSIDES_FLEXIBLEBILFORMRHS_H_
#define LAWA_METHODS_RB_RIGHTHANDSIDES_FLEXIBLEBILFORMRHS_H_

#include <vector>
#include <lawa/methods/adaptive/datastructures/datastructures.h>


namespace lawa {

template <typename Index, typename LocalOperatorType>
class FlexibleBilformRhs {

public:

	typedef typename LocalOperatorType::T T;

    FlexibleBilformRhs(){};

	FlexibleBilformRhs(std::vector<LocalOperatorType*>& _bilformvec);

	virtual
	T
    operator()(const Index &index);

	virtual
	Coefficients<Lexicographical,T,Index>
	operator()(const IndexSet<Index> &indexset);

	virtual
	void
	set_active_comp(int i);

	virtual
	void
	set_active_u(Coefficients<Lexicographical,T,Index> const* _u);

    const std::vector<LocalOperatorType*>&
    get_bilformvec() const;
	
	void
    clear();

protected:

    std::vector<LocalOperatorType*>& 				bilformvec;
    std::vector<int> 								active_comp;

    Coefficients<Lexicographical,T,Index> const* 	active_u;

	FlexibleBilformRhs(const FlexibleBilformRhs& rhs);
};

} // namespace lawa

#include <lawa/methods/rb/righthandsides/flexiblebilformrhs.tcc>

#endif /* LAWA_METHODS_RB_RIGHTHANDSIDES_FLEXIBLEBILFORMRHS_H_ */
