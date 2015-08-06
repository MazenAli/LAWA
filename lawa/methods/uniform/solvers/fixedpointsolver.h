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
 
#ifndef LAWA_METHODS_UNIFORM_SOLVERS_FIXEDPOINTSOLVER_H
#define LAWA_METHODS_UNIFORM_SOLVERS_FIXEDPOINTSOLVER_H 1

#include <lawa/flensforlawa.h>

namespace lawa{
    
/* Fixed Point Solver for periodic problems
 *      This class solves periodic problems by repeatedly calling
 *      an underlying method (e.g. TimeStepping) with u_0^(i) = u(T)^(i-1)
 *      until || u_0^(i) - u(T)^(i) ||_l2 < tol 
 */
template<typename T, typename Method>
class FixedPointSolver
{
    public:
        typedef typename Method::RHSType RHSType;
        
        FixedPointSolver(Method& _method);
        
        flens::DenseVector<flens::Array<T> >
        solve(flens::DenseVector<flens::Array<T> > u_0, bool saveSols = false, 
              int maxIterations = 1000, T tol = 10e-15);
        
        flens::DenseVector<flens::Array<T> >
        solve(flens::DenseVector<flens::Array<T> > u_0,
                flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& fmatrix,
                int maxIterations = 1000, T tol = 10e-15);
        
        flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >&
        getSolutions(){ return method.getSolutions();}
    
        void setLevel(int level){ method.setLevel(level);}
        void setRHS(RHSType& rhs);
        
        
    private:
        Method& method;
        
        T
        getError(flens::DenseVector<flens::Array<T> >& u1, flens::DenseVector<flens::Array<T> >& u2);
    
};    
  
} // namespace lawa

#include <lawa/methods/uniform/solvers/fixedpointsolver.tcc>

#endif // LAWA_METHODS_UNIFORM_SOLVERS_FIXEDPOINTSOLVER_H

