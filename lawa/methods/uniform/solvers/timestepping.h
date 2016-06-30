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
 
#ifndef LAWA_METHODS_UNIFORM_SOLVERS_TIMESTEPPING_H
#define LAWA_METHODS_UNIFORM_SOLVERS_TIMESTEPPING_H 1

namespace lawa{

/* TimeStepping
 *      This class performs a time-stepping method.
 *      For each t=1,..,K, the underlying solver (e.g. a ThetaScheme) is called.
 *      Optionally, the solutions for each t are stored in U (u(t_k) in k-th col).
 */
template <typename T, typename Solver>
class TimeStepping
{
    public:
        
        typedef typename Solver::RHSType RHSType;
        
        TimeStepping(Solver& _solver, T _deltaT, FLENS_DEFAULT_INDEXTYPE _timesteps, FLENS_DEFAULT_INDEXTYPE _levelX);

        flens::DenseVector<flens::Array<T> > 
        solve(flens::DenseVector<flens::Array<T> >& u_0, bool saveSols = false);
    
        flens::DenseVector<flens::Array<T> > 
        solve(flens::DenseVector<flens::Array<T> >& u_0, 
              flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& fmatrix);
    
        flens::DenseVector<flens::Array<T> > 
        getResiduum(flens::DenseVector<flens::Array<T> >& u);

        T getDeltaT(){ return deltaT;}
        T getSteps(){ return timesteps;}
        T getLevel(){ return levelX;}
        
        flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >&
        getSolutions();
               
        flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >&
        getSolutions(flens::DenseVector<flens::Array<T> >& u);
           
        void setDeltaT(T delta){ deltaT = delta;}
        void setSteps(FLENS_DEFAULT_INDEXTYPE steps){ timesteps = steps;}
        void setLevel(FLENS_DEFAULT_INDEXTYPE J){ levelX = J;}
        
        void setRHS(RHSType& rhs){ solver.setRHS(rhs);}
        
    private:
        Solver& solver;
        T deltaT;
        FLENS_DEFAULT_INDEXTYPE timesteps;
        FLENS_DEFAULT_INDEXTYPE levelX;
        
        flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > U;

};

} // namespace lawa

#include <lawa/methods/uniform/solvers/timestepping.tcc>

#endif // LAWA_METHODS_UNIFORM_SOLVERS_TIMESTEPPING_H

