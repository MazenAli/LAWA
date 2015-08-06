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

#ifndef LAWA_METHODS_RB_DATASTRUCTURES_THETASTRUCTURE_H_
#define LAWA_METHODS_RB_DATASTRUCTURES_THETASTRUCTURE_H_

#include <vector>
#include <array>

namespace lawa {

/* ThetaStructure:
 * 		Manage a vector of theta functions.
 * 		Helpful in order to store them in only one place.
 *
 * 		Addionally a pointer to the current parameter, so that
 * 		affine operators/rhs that only have an operator()(Index)
 * 		can still evaluate the functions
 */
template<typename ParamType>
class ThetaStructure{

public:
	typedef typename ParamType::value_type T;
    typedef T (*ThetaFct)(const ParamType& param); // Argumente -> eher auch RBThetaData-Objekt?

    ThetaStructure();

    ThetaStructure(const std::vector<ThetaFct>& _thetas);

    size_t
    size() const;

    void
    set_param(const ParamType& _param);

    ParamType&
    get_param();

    T
    eval(size_t i, const ParamType& mu) const;

    T
    eval(size_t i) const;

private:
	std::vector<ThetaFct> 	thetas;
	ParamType    			current_param;

	ThetaStructure(const ThetaStructure& thetastructure);
};

} // namespace lawa

#include <lawa/methods/rb/datastructures/thetastructure.tcc>

#endif /* LAWA_METHODS_RB_DATASTRUCTURES_THETASTRUCTURE_H_ */
