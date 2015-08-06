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
 
#ifndef LAWA_METHODS_UNIFORM_SOLVERS_MULTIGRIDSOLVER_H
#define LAWA_METHODS_UNIFORM_SOLVERS_MULTIGRIDSOLVER_H 1

namespace lawa{

/* Multigrid Solver
 *      This class provides the framework for multigrid methods by
 *      implementing a general V- and W-Cycle for Wavelets 
 *      (using decompose/reconstruct for projections and restrictions) 
 */
template<typename T, typename PrimalBasis, typename DualBasis, typename Smoother, typename Solver>
class MultigridSolver{
    
      typedef flens::DenseVector<flens::Array<T> > DenseVectorT;
    
  private:
      PrimalBasis& primalbasis;
      DualBasis& dualbasis;
      Smoother& smoother;
      Solver& solver; 
      int nu1, nu2;                 
      int minLevel;
    
  public:
      MultigridSolver(PrimalBasis& _primalbasis, DualBasis& _dualbasis, Smoother& _smoother, 
                      Solver& _solver, int _nu1, int _nu2, int _minLevel = 0);
      
      DenseVectorT
      vCycle(int i, int level, DenseVectorT& u, DenseVectorT& f);
                  
      DenseVectorT
      wCycle(int i, int level, DenseVectorT& u, DenseVectorT& f);      
      
      int getMinLevel(){ return minLevel;}
};

} //  namespace lawa

#include <lawa/methods/uniform/solvers/multigridsolver.tcc>

#endif // LAWA_METHODS_UNIFORM_SOLVERS_MULTIGRIDSOLVER_H

