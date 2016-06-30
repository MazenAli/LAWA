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

#ifndef LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_FLEXIBLECOMPOUNDRHS_H
#define LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_FLEXIBLECOMPOUNDRHS_H 1

#include <vector>
#include <lawa/methods/adaptive/datastructures/datastructures.h>

namespace lawa {

template <typename T, typename Index, typename RHSType>
class FlexibleCompoundRhs {

public:

	FlexibleCompoundRhs(std::vector<RHSType*>& _rhsvec);

	virtual
	T
    operator()(const Index &index);

	virtual
	Coefficients<Lexicographical,T,Index>
	operator()(const IndexSet<Index> &indexset);

    template <typename Prec>
	Coefficients<Lexicographical,T,Index>
	operator()(const IndexSet<Index> &indexset, Prec& P);
	
    virtual
	void
	set_active_comp(FLENS_DEFAULT_INDEXTYPE i);

    virtual
    const RHSType&
    get_comp(FLENS_DEFAULT_INDEXTYPE i) const;
	
	void
    clear();

protected:

    std::vector<RHSType*>& 	rhsvec;
    std::vector<FLENS_DEFAULT_INDEXTYPE> 		active_comp;

private:

    FlexibleCompoundRhs(const FlexibleCompoundRhs& rhs);
};

}   // namespace lawa

#include <lawa/methods/adaptive/righthandsides/flexiblecompoundrhs.tcc>

#endif  // LAWA_METHODS_ADAPTIVE_RIGHTHANDSIDES_FLEXIBLECOMPOUNDRHS_H
