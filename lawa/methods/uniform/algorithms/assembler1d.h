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
 
#ifndef LAWA_METHODS_UNIFORM_ALGORITHMS_ASSEMBLER1D_H
#define LAWA_METHODS_UNIFORM_ALGORITHMS_ASSEMBLER1D_H 1

#include <extensions/extensions.h>

namespace lawa{    

template<typename T, typename Basis>
class Assembler1D
{
    private:
        const Basis& basis;
   
    public: 
        Assembler1D(const Basis& _basis);

        /* Assemble Stiffness Matrix in transposed form, i.e.
         * corresponding to a(v,u)
         */
        template <typename BilinearForm>
        flens::SparseGeMatrix<flens::extensions::CRS<T,flens::CRS_General> >
        assembleStiffnessMatrix(BilinearForm& a, FLENS_DEFAULT_INDEXTYPE J, T tol = 10e-15);
    
        template <typename RHSIntegral>
        flens::DenseVector<flens::Array<T> >
        assembleRHS(RHSIntegral& rhs, FLENS_DEFAULT_INDEXTYPE J);
        
        template <typename Preconditioner>
        flens::DiagonalMatrix<T>    
        assemblePreconditioner(Preconditioner& P, FLENS_DEFAULT_INDEXTYPE J);
};

} // namespace lawa

#include <lawa/methods/uniform/algorithms/assembler1d.tcc>

#endif // LAWA_METHODS_UNIFORM_ALGORITHMS_ASSEMBLER1D_H

