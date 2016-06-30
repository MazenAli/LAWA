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


#ifndef LAWA_OPERATORS_PDEOPERATORS1D_WEIGHTEDCONVECTIONOPERATOR1D_H
#define LAWA_OPERATORS_PDEOPERATORS1D_WEIGHTEDCONVECTIONOPERATOR1D_H 1

#include <lawa/settings/enum.h>
#include <lawa/functiontypes/function.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/integrals/integral.h>

namespace lawa {

/* Weighted Convection OPERATOR 1D
 *
 *    a(v,u) =  Integral(w(x) * v * u_x)
 *
 */
template <typename T, typename Basis, QuadratureType Quad>
class WeightedConvectionOperator1D{
    
    public:
        
        const Basis& basis;

        WeightedConvectionOperator1D(const Basis& _basis, Function<T> weightFct);

        T
        operator()(XType xtype1, FLENS_DEFAULT_INDEXTYPE j1, FLENS_DEFAULT_INDEXTYPE k1,
                   XType xtype2, FLENS_DEFAULT_INDEXTYPE j2, FLENS_DEFAULT_INDEXTYPE k2) const;

        T
        operator()(const Index1D &row_index, const Index1D &col_index) const;
        
    private:
    
        Function<T> W;

        IntegralF<Quad, Basis, Basis> integral;

};

} // namespace lawa

#include <lawa/operators/pdeoperators1d/weightedconvectionoperator1d.tcc>

#endif  // LAWA_OPERATORS_PDEOPERATORS1D_WEIGHTEDCONVECTIONOPERATOR1D_H

