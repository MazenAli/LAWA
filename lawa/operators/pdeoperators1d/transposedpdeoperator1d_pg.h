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

#ifndef LAWA_OPERATORS_PDEOPERATORS1D_TRANSPOSEDPDEOPERATOR1D_PG_H
#define LAWA_OPERATORS_PDEOPERATORS1D_TRANSPOSEDPDEOPERATOR1D_PG_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/integrals/integral.h>
#include <lawa/settings/enum.h>

namespace lawa {

/* Transposed PDE ConstCoefficient OPERATOR : Petrov-Galerkin version
 *
 *    a^T(v,u) =  [diffusion * Integral(v_x * u_x) +  convection * Integral(v * u_x)
 *              + reaction * Integral(v * u)]^T
 *
 */
template <typename T, typename TrialBasis, typename TestBasis>
class TransposedPDEOperator1D_PG{
    
    public:

    	const TrialBasis& trialbasis;
    	const TestBasis& testbasis;
        T reaction, convection, diffusion;

        TransposedPDEOperator1D_PG(const TrialBasis& _trialbasis, const TestBasis& _testbasis, T _reaction, T _convection, T _diffusion);

        T
        operator()(XType xtype1, int j1, long k1,
                   XType xtype2, int j2, long k2) const;

        T
        operator()(const Index1D &row_index, const Index1D &col_index) const;
        
    private:

        Integral<Gauss, TestBasis, TrialBasis> integral;
    

};

}   //namespace lawa

#include <lawa/operators/pdeoperators1d/transposedpdeoperator1d_pg.tcc>

#endif  // LAWA_OPERATORS_PDEOPERATORS1D_TRANSPOSEDPDEOPERATOR1D_PG_H

