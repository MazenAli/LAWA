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

#ifndef LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_ABSTRACTLOCALOPERATOR1D_H_
#define LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_ABSTRACTLOCALOPERATOR1D_H_

#include <lawa/methods/adaptive/datastructures/treecoefficients1d.h>

namespace lawa {

template <typename _T, typename TestBasis, typename TrialBasis=TestBasis>
class AbstractLocalOperator1D {
public:
    typedef _T T;

    virtual
    void
    eval(const TreeCoefficients1D<T> &PsiLambdaHat,
         TreeCoefficients1D<T> &PsiLambdaCheck,
         const char* mode) = 0;

    virtual
    const TrialBasis&
    getTrialBasis() const  = 0;

    virtual
    const TestBasis&
    getTestBasis() const   = 0;
};


} // namespace lawa

#endif /* LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_ABSTRACTLOCALOPERATOR1D_H_ */
